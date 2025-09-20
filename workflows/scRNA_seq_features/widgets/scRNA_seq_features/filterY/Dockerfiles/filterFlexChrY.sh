#!/bin/bash
function  checkPids() {
    for pid in "${pidList[@]}"; do
        if ! kill -0 "$pid" 2>/dev/null; then
            # If PID is no longer running, remove it from the list
            pidList=(${pidList[@]/$pid})
            return 0
        fi
    done
    return 1
}
function filterBam() {
    # Split a BAM into reads that ARE and are NOT present in a read-name list
    # using a single, multi-threaded samtools pass.
    #
    # Usage: filterBam <input.bam> <read_names.txt> <output_Y.bam> <output_noY.bam> [threads]

    local input_bam="$1"
    local ynames_file="$2"
    local output_bam_Y="$3"
    local output_bam_noY="$4"
    local threads="${5:-4}"  # default to 4 compression threads

    if [[ $# -lt 4 || $# -gt 5 ]]; then
        echo "Usage: filterBam <input.bam> <read_names.txt> <output_Y.bam> <output_noY.bam> [threads]" >&2
        return 1
    fi

    # Ensure samtools is available
    if ! command -v samtools &>/dev/null; then
        echo "Error: 'samtools' not found in PATH." >&2
        return 1
    fi

    # Validate input files
    [[ -f "$input_bam" ]] || { echo "Error: Input BAM '$input_bam' not found." >&2; return 1; }
    [[ -f "$ynames_file" ]] || { echo "Error: Read-name list '$ynames_file' not found." >&2; return 1; }

    echo "→ Filtering '$input_bam' with $(wc -l < "$ynames_file") read names (threads=$threads)…"

    samtools view -@ "$threads" \
                  -N "$ynames_file" \
                  -U "$output_bam_noY" \
                  -b "$input_bam" > "$output_bam_Y"

    local status=$?
    if [[ $status -ne 0 ]]; then
        echo "samtools view failed with status $status" >&2
        return $status
    fi

    echo "   Y-matching reads written to:     $output_bam_Y"
    echo "   non-Y reads written to:          $output_bam_noY"
}
function filter_fastqs_in_directory() {
    local fastq_dir="$1"
    local ynames_file="$2"
    local output_dir="$3"
    
    # Collect all FASTQ files in the directory
    mapfile -t fastqs < <(
        find "$fastq_dir" -type f \
            \( -name "*.fq" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fastq.gz" \)
    )
    
    echo "→ Found ${#fastqs[@]} FASTQ files in $fastq_dir – processing with $nThreads parallel jobs"
    
    # Launch parallel background jobs with throttling
    local pidList=()
    for fastq in "${fastqs[@]}"; do
        (
            # Create output file names in subshell
            local fq_name
            fq_name=$(basename "$fastq")
            
            # Remove extension to build output names
            local base_name="${fq_name%.*}"
            # Handle .fastq.gz -> remove both extensions
            [[ "$fq_name" == *.fastq.gz ]] && base_name="${fq_name%.fastq.gz}"
            # Handle .fq.gz -> remove both extensions  
            [[ "$fq_name" == *.fq.gz ]] && base_name="${fq_name%.fq.gz}"
            
            local output_fastq_Y="${output_dir}/${base_name}_Y.fastq"
            local output_fastq_noY="${output_dir}/${base_name}_noY.fastq"
            
            filter_fastq "$fastq" "$ynames_file" "$output_fastq_Y" "$output_fastq_noY"
        ) &
        
        pidList+=($!)
        
        # Throttle to nThreads concurrent jobs
        while (( ${#pidList[@]} >= nThreads )); do
            checkPids || sleep 0.5
        done
    done
    
    # Wait for all remaining jobs to complete
    for pid in "${pidList[@]}"; do
        wait "$pid"
    done
    
    echo "✓ All FASTQ filtering complete"
}
function filter_fastq() {
    # --- 1. Argument Validation ---
    local input_fastq="$1"
    local read_list="$2"
    local out_in_list="$3"
    local out_not_in_list="$4"

    if [[ $# -ne 4 ]]; then
        echo "Error: Invalid number of arguments." >&2
        echo "Usage: filter_fastq <input.fastq[.gz]> <read_names.txt> <output_in.fastq> <output_not_in.fastq>" >&2
        return 1
    fi

    if [[ ! -f "$input_fastq" ]]; then
        echo "Error: Input FASTQ file not found at '$input_fastq'" >&2
        return 1
    fi

    if [[ ! -f "$read_list" ]]; then
        echo "Error: Read list file not found at '$read_list'" >&2
        return 1
    fi

    # --- 2. Determine Compression and Set Up Commands ---
    local output_suffix=""
    local input_reader_cmd="cat"
    local output_writer_cmd="cat"

    # Check if the input file ends with .gz
    if [[ "${input_fastq##*.}" == "gz" ]]; then
        echo "Compressed input detected. Output will be gzipped."
        output_suffix=".gz"
        input_reader_cmd="gunzip -c"
        output_writer_cmd="gzip"
    else
        echo "Uncompressed input detected."
    fi

    # Define the full output filenames and the commands to write to them
    local final_out_in="${out_in_list}${output_suffix}"
    local final_out_not_in="${out_not_in_list}${output_suffix}"
    local write_in_cmd="${output_writer_cmd} > '${final_out_in}'"
    local write_not_in_cmd="${output_writer_cmd} > '${final_out_not_in}'"

    # --- 3. Define the Awk Script ---
    # This script is passed as a single string to the awk command.
    # It uses several awk features for efficiency:
    #
    # - FNR==NR: This pattern is only true while awk is reading the first file
    #   it's given (the read_list). It loads every line from that file into an
    #   associative array called `reads` for fast lookups later.
    #
    # - FNR % 4 == 1: In a FASTQ file, the read header is on every 4th line
    #   starting with line 1. When this pattern matches, we know we are on a
    #   header line. The script then checks if this read name exists in our
    #   `reads` array and sets a flag `in_list`.
    #
    # - The main block `{...}` runs for every line. It prints the current line
    #   to one of two pipes ('cmd_in' or 'cmd_not_in') based on the `in_list` flag.
    #
    # - The END block ensures the file handles for the output pipes are closed
    #   cleanly, which is crucial for piped commands like gzip.
    #
    local awk_script='
        BEGIN {
            # Receive the shell command strings into awk variables
            cmd_in = "'"${write_in_cmd}"'"
            cmd_not_in = "'"${write_not_in_cmd}"'"
        }

        # Part 1: Process the read_list file
        FNR==NR {
            reads[$1]=1
            next
        }

        # Part 2: Process the FASTQ file
        # This is the header line for a read
        FNR % 4 == 1 {
            # Extract readname from the first field (e.g., @SEQ_ID/1)
            readname = $1
            # Remove leading "@"
            sub(/^@/, "", readname)
            # Remove trailing pair info (e.g., /1, /2, or " 1", " 2")
            sub(/[ \t\/].*$/, "", readname)
            # Set a flag indicating if this read should be kept
            in_list = (readname in reads)
        }

        # This block executes for all 4 lines of a read
        {
            if (in_list) {
                print | cmd_in
            } else {
                print | cmd_not_in
            }
        }

        END {
            close(cmd_in)
            close(cmd_not_in)
        }
    '

    # --- 4. Execute the Process ---
    echo "Starting FASTQ filtering..."
    echo "  - Reads in list will be written to: ${final_out_in}"
    echo "  - Reads not in list will be written to: ${final_out_not_in}"

    # The input reader command (cat or gunzip -c) streams the FASTQ file to stdout.
    # This is piped into awk.
    # Awk is given the script, the read_list file, and '-' which tells it to
    # read the FASTQ data from stdin (the pipe).
    $input_reader_cmd "$input_fastq" | awk "$awk_script" "$read_list" -

    echo "Filtering complete."
}

# Set default number of parallel threads for FASTQ processing
: "${nThreads:=16}"

alignedDir=/storage/scRNAseq_output/Alignments
sequenceDir=/storage/JAX_scRNAseq01
outputDir=/storage/Y_filtered

[ -z "$alignedDir" ] && echo "need $alignedDir" && exit 1
[ -z "$sequenceDir" ] && echo "need $sequenceDir" && exit 1
[ -z "$outputDir" ] &&  echo "need $outputDir" && exit 1
[ -z "$aligner" ] && aligner=star


trimmedNames=($(find $alignedDir -mindepth 1 -maxdepth 1 -type d))

[ -n "$sequenceDir" ] &&  sequenceFiles=($(find $sequenceDir -type f \( -name "*.fq" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fastq.gz" \)))

for trimmedName in "${trimmedNames[@]}"; do
    echo "working on $trimmedName"
    output="$outputDir/$(basename $trimmedName)/$aligner"
    mkdir -p "$output"
    ynames="$output/ynames.txt"
    bamIn="$trimmedName/$aligner/Aligned.sortedByCoord.out.bam"
    bamYOut="$output/Aligned.sortedByCoord_Y_.out.bam"
    bamNoYOut="$output/Aligned.sortedByCoord_noY_.out.bam"
    #find readnames that align to Y chromosome    
    [[ ! -f "$ynames" || -n "$overwrite" ]] && \
        echo "samtools view "$bamIn" | awk '\$3 == \"chrY\" { unique[\$1]++ } END { for (val in unique) print val }' > $ynames" && \
        eval "samtools view "$bamIn" | awk '\$3 == \"chrY\" { unique[\$1]++ } END { for (val in unique) print val }' > $ynames"
    #check if the ynames file is empty
    if [ ! -s "$ynames" ]; then
        echo "No Y chromosome reads found in $bamIn"
        continue
    fi
    echo "filtering bam files"
    #filter the bam files unless the output files already exist and overwrite is set to true
    if [ ! -f "$bamYOut" ] || [ ! -f "$bamNoYOut" ] || [ -n "$overwrite" ]; then
        echo "filterBam $bamIn $ynames $bamYOut $bamNoYOut"
        filterBam $bamIn $ynames $bamYOut $bamNoYOut
    fi

    #filter the fastqs
    filter_fastqs_in_directory $sequenceDir $ynames $output
done
