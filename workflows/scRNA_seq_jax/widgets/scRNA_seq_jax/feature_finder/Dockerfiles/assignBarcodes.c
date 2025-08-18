#include <ctype.h>
#include <getopt.h>  // For command line parsing
#include <glib.h>
#include <limits.h>  // For CHAR_BIT
#include <math.h>
#include <omp.h>
#include <stddef.h>  // For offsetof
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> // For mkdir and file exists
#include <sys/time.h> //for timing
#include <unistd.h>
#include <zlib.h>
#include "defines.h"
#include "queue.h"

//code for feature sequences stats

typedef struct _feature_arrays {
    int number_of_features;
    int max_length;
    char **feature_names;
    char *feature_names_storage;
    unsigned int *feature_lengths;
    unsigned char *feature_code_lengths;
    char **feature_sequences;
    char *feature_sequences_storage;
    unsigned char **feature_codes;
    unsigned char *feature_codes_storage; 
} feature_arrays;

typedef struct _feature_counts {
    unsigned char sequence_code[4]; // need to hardcode to 4 for 32 bit hash - if BARCODE_CODE_LENGTH > 4 then need to change this
    uint16_t counts[];
} feature_counts;

typedef struct _feature_sequences {
    uint32_t counts;
    char hamming_distance;
    unsigned char feature_index;
    char sequence[];
} feature_sequences;

typedef struct _feature_umi_counts {
    unsigned char sequence_umi_code[8]; // hard coded to 8 for 64 bit
    uint16_t counts[];
} feature_umi_counts;

typedef struct _memory_pool {
    size_t block_size;       // Size of each block (size of features_block)
    size_t blocks_per_pool;  // Number of blocks per pool
    size_t free_blocks;      // Number of free blocks
    void *next_free;         // Pointer to the next free block
    struct _storage_block *first_block;   // Pointer to the first memory block
    struct _storage_block *current_block; // Pointer to the current memory block
} memory_pool;

typedef struct _memory_pool_collection {
    memory_pool *feature_counts_pool;
    memory_pool *feature_umi_counts_pool;
    memory_pool *feature_sequences_pool;
    memory_pool *unmatched_barcodes_features_block_pool;
} memory_pool_collection;

typedef struct _sample {
    int nFiles;
    char **fastqFileR1;
    char **fastqFileR2;
} sample;

typedef struct _storage_block {
    struct _storage_block *next;
    unsigned char *storage;
} storage_block;

typedef struct _unmatched_barcodes_features {
    struct _unmatched_barcodes_features_block *next;
    unsigned char feature_index;
    unsigned char *barcode;
    unsigned char *umi;
    unsigned char number_of_closest_barcodes;
    unsigned char *closest_barcodes;
    unsigned char *Qscores;
} unmatched_barcodes_features;

typedef struct _unmatched_barcodes_features_block {
    struct _unmatched_barcodes_features_block *next;
    unsigned char storage[];
} unmatched_barcodes_features_block;

typedef struct _unmatched_barcodes_features_block_list {
    unmatched_barcodes_features_block *first_entry;
    unmatched_barcodes_features_block *last_entry;
} unmatched_barcodes_features_block_list;

typedef struct _unit_sizes {
    // stores unit sizes for dynamically allocated structs
    size_t feature_counts;
    size_t feature_umi_counts;
    size_t feature_sequences;
    size_t unmatched_barcodes_features_block;
} unit_sizes;

typedef struct _statistics {//keep track of stats
    double start_time;
    size_t nMismatches;
    size_t recovered;
    size_t pending;
    size_t valid;
    size_t pending_recovered;  
    size_t total_unmatched_features;
    size_t number_of_reads;
    unmatched_barcodes_features_block_list unmatched_list;
} statistics;

typedef struct _data_structures {
    GHashTable *filtered_hash;
    GHashTable *sequence_umi_hash;
    GHashTable *unique_features_match;
    Queue *neighbors_queue;
} data_structures;

//function prototypes sorted alphabetically
void add_deduped_count(feature_counts *s, uint32_t *counts, uint16_t stringency, uint32_t min_counts);
unmatched_barcodes_features_block* add_unmatched_barcode_store_feature(unsigned char *barcodes, char *umi, unsigned char *qscores, unsigned char feature_index, int number_of_variants, memory_pool_collection *pools, statistics *stats);
storage_block* allocate_storage_block(size_t block_size);
char check_neighbor(uint64_t code64, uint32_t *counts, data_structures *hashes);
int checkAndCorrectBarcode(char **lines, int maxN, unsigned char feature_index, data_structures *hashes, memory_pool_collection *pools, statistics *stats);
int checkAndCorrectFeature(char *line, feature_arrays *features, int maxHammingDistance, int nThreads, int *hamming_distance, char *matching_sequence, int maxN, char *ambiguous);
char check_sequence(char *sequence, int sequence_length);
void code2string(unsigned char *code, char *string, int length);
int compare_feature_sequences(const void *a, const void *b);
void destroy_data_structures(data_structures *hashes);
void initialize_data_structures(data_structures *hashes);
int existing_output_skip(char keep_existing, char *directory);
void expand_memory_pool(memory_pool *pool);
void finalize_processing(feature_arrays *features, data_structures *hashes,  char *directory, int debug, memory_pool_collection *pools, statistics *stats);
unsigned char* find_best_posterior_match (unmatched_barcodes_features_block *entry_block, int number_of_features, double min_posterior,statistics *stats);
int find_closest_barcodes(unsigned char* code, unsigned char *corrected_codes, unsigned char *indices);
void find_connected_component(gpointer start_key, uint32_t *counts, data_structures *hashes);
void find_deduped_counts(data_structures *hashes);
int find_feature_match_parallel(feature_arrays *features, char *lineR2, int maxHammingDistance, int nThreads, int *bestScore, char **matching_sequence);
int find_feature_match_single(feature_arrays *features, char *lineR2, int maxHammingDistance,int *bestScore, char **matching_sequence);
int find_neighbors(uint64_t key64, uint64_t *neighbors, uint32_t *counts, data_structures *hashes);
int find_matching_codes(unsigned char *sequence_code, unsigned char *feature_code, size_t length, size_t feature_length);
void find_sample_directory_name(char *directory, sample *sample, char *sample_directory, char sample_flag);
int find_substring_match(const char *query, feature_arrays *features);
int find_variant_match(unsigned char *code, int sequence_index, unsigned char *corrected_codes);
double get_time_in_seconds();
int generate_sequences(char* string, int string_length, int *indices, char **output, int index, int number_of_indices, char *output_indices[], int number_of_outputs);
void initcode2seq();
void initdiff2hamming(unsigned char *difference);
void initseq2Code();
void initialize_statistics(statistics *stats);
memory_pool* initialize_storage_pool(size_t block_size, size_t blocks_per_pool);
void initialize_unit_sizes();
int insert_feature_sequence(char *sequence, unsigned char feature_index, unsigned char hamming_distance, data_structures *hashes, memory_pool_collection *pools);
int mkdir_p(const char *path);
void open_R1R2_files(const char *fastqFileR1, const char *fastqFileR2, gzFile *fastqFileR1gz, gzFile *fastqFileR2gz);
void printFeatureCounts(feature_arrays *custom_features, int *deduped_counts, int *barcoded_counts, int **coexpression_counts, int **coexpression_histograms, char *directory, data_structures *hashes, statistics *stats);
int print_feature_sequences(feature_arrays *features, int *total_counts, char *directory, data_structures *hashes);
void process_pending_barcodes(data_structures *hashes, memory_pool_collection *pools, statistics *stats);
void process_reads( data_structures *hashes, char *lines[8], feature_arrays *features, int maxHammingDistance, int nThreads, memory_pool_collection *pools, statistics *stats);
void read_unmatched_features_block(unmatched_barcodes_features_block *entry_block, unmatched_barcodes_features *entry);
unsigned char* readWhiteList(char *whitelist_filename, GHashTable *hash);
feature_arrays* readTagsFile(const char* filename);
char* printToBuffer(char *string, int sequence_length, char *buffer);
int string2code(char *string, int sequence_length, unsigned char *code);
int string2code_debug(char *string, int sequence_length, unsigned char *code);
void string2all_codes(char *string, unsigned char codes[][LINE_LENGTH/2+1], int *lengths);
void update_feature_counts(char *barcode_string, char *umi, uint32_t feature_index, data_structures *hashes,memory_pool_collection *pools);
void update_feature_counts_from_code(unsigned char *code, char *umi, uint32_t feature_index, data_structures *hashes, memory_pool_collection *pools);
void update_umi_counts(unsigned char *code, char *umi, uint32_t feature_index, data_structures *hashes, memory_pool_collection *pools);
void process_files_in_sample(char *directory, char **fastqFileR1, char **fastqFileR2, int nR1files, feature_arrays *features, int maxHammingDistance, int nThreads, int debug, memory_pool_collection *pools, statistics *stats, data_structures *hashes);
void get_fastq_filenames(int positional_arg_count, int argc, char *fastqFileR2[], char *argv[], char *fastqFileR1[], char *fastqFilesR1string, char *fastqFilesR2string);
int organize_fastq_files_into_samples( char *fastqFileR1[], char *fastqFileR2[], int nR1files, sample **samples);
void add_sample(sample *entry, char *fastqFileR1, char *fastqFileR2);
void free_memory_pool_collection(memory_pool_collection *pools);
void free_memory_pool_storage(memory_pool *pool);

//tables for converting sequences to codes
static unsigned char seq2code[256];
static char code2seq[256][4];
static unsigned char diff2Hamming[256];

//unit sizes for dynamically allocated structs
static unit_sizes dynamic_struct_sizes;



//global variables for the program

static int debug;

//valid barcode list and hash
unsigned char *whitelist;
GHashTable *whitelist_hash; 
//hash tables for storing feature counts, umi counts, feature sequences, whitelist and filtered barcodes

//size globals to replace constants
static int barcode_length;
static int barcode_code_length;
static int number_of_features;
static int maximum_feature_length;
static int feature_code_length;

//default values for the program
static int max_feature_n=MAX_FEATURE_N;
static int max_barcode_n=MAX_BARCODE_N;
static int max_barcode_mismatches=MAX_BARCODE_MISMATCHES;
static int umi_length=UMI_LENGTH;
static int umi_code_length=UMI_CODE_LENGTH;
static int constant_offset=-1;
static uint16_t stringency=500;
static uint32_t min_counts=0;
void initialize_statistics(statistics *stats) {
    stats->start_time = get_time_in_seconds();
    stats->nMismatches = 0;
    stats->recovered = 0;
    stats->pending = 0;
    stats->valid = 0;
    stats->pending_recovered = 0;
    stats->total_unmatched_features = 0;
    stats->number_of_reads = 0;
    stats->unmatched_list.first_entry = NULL;
    stats->unmatched_list.last_entry = NULL;
}


void free_memory_pool_storage(memory_pool *pool) {
    storage_block *current = pool->first_block;
    while (current != NULL) {
        storage_block *next = current->next;
        free(current->storage);
        free(current);
        current = next;
    }
}

void free_unmatched_barcodes_features_list(unmatched_barcodes_features_block_list *list) {
    unmatched_barcodes_features_block *current = list->first_entry;
    while (current != NULL) {
        unmatched_barcodes_features_block *next = current->next;
        free(current);
        current = next;
    }
    list->first_entry = NULL;
    list->last_entry = NULL;
}


double get_time_in_seconds() {
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec + (time.tv_usec / 1000000.0);
}
int mkdir_p(const char *path) {
    char temp[1024];
    char *p = NULL;
    size_t len;

    // Copy path and ensure it ends with '/'
    snprintf(temp, sizeof(temp), "%s", path);
    len = strlen(temp);
    if (temp[len - 1] == '/')
        temp[len - 1] = 0;

    // Iterate through each directory in the path
    for (p = temp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;

            // Create directory if it doesn't exist
            if (mkdir(temp, S_IRWXU) != 0 && errno != EEXIST) {
                perror("mkdir");
                return -1;
            }
            *p = '/';
        }
    }
    // Create the final directory
    if (mkdir(temp, S_IRWXU) != 0 && errno != EEXIST) {
        perror("mkdir");
        return -1;
    }
    return 0;
}

void lowerCaseDifferences(char *ref, char *query, int length){
    for (int i=0; i<length; i++){
        if (ref[i] != query[i]){
            query[i]=tolower(query[i]);
        }
    } 
}   
int compare_feature_sequences(const void *a, const void *b) {
    const feature_sequences *fa = *(const feature_sequences **)a;
    const feature_sequences *fb = *(const feature_sequences **)b;

    // First sort by feature_index
    if (fa->feature_index != fb->feature_index)
        return fa->feature_index - fb->feature_index;

    // Then sort by hamming_distance
    if (fa->hamming_distance != fb->hamming_distance)
        return fa->hamming_distance - fb->hamming_distance;

    // Lastly sort by counts
    return fb->counts - fa->counts;
}

void initseq2Code(){
    for (int i=0; i<256; i++){
        seq2code[i]=0;
    }
    seq2code['A']=0;
    seq2code['C']=1;
    seq2code['G']=2;
    seq2code['T']=3;
}
void initcode2seq(){
    char temp[1025];
    strcpy(temp, LOOKUP_STRING);
    for (int i=0; i<256; i++){
        code2seq[i][0]=temp[4*i];
        code2seq[i][1]=temp[4*i+1];
        code2seq[i][2]=temp[4*i+2];
        code2seq[i][3]=temp[4*i+3];
    }
}
void initdiff2hamming(unsigned char *difference){
    memset(difference,0,256);
    for (int i=0; i<256; i++){
        difference[i]=0;
          unsigned char mask=3;
        //check if first 2 bits are 0
        for (int j=0; j<4; j++){
            if (i & mask){
                difference[i]++;
            }
            mask=mask << 2;
        }
        //DEBUG_PRINT( "Difference %d %d\n", i, difference[i]);
    }
}
memory_pool_collection* initialize_memory_pool_collection() {
    memory_pool_collection *pools = (memory_pool_collection *)malloc(sizeof(memory_pool_collection));
    if (!pools) {
        perror("Failed to allocate memory for memory pool collection");
        exit(EXIT_FAILURE);
    }
    pools->feature_counts_pool=initialize_storage_pool(dynamic_struct_sizes.feature_counts, FEATURE_BARCODE_BLOCK_SIZE);
    pools->feature_umi_counts_pool=initialize_storage_pool(dynamic_struct_sizes.feature_umi_counts, UMI_STORAGE_BLOCK);
    pools->feature_sequences_pool=initialize_storage_pool( dynamic_struct_sizes.feature_sequences, FEATURE_SEQUENCE_BLOCK_SIZE);
    pools->unmatched_barcodes_features_block_pool=initialize_storage_pool(dynamic_struct_sizes.unmatched_barcodes_features_block, BARCODE_STORAGE_BLOCK_SIZE);
    return pools;
}

void free_memory_pool_collection(memory_pool_collection *pools) {
    free_memory_pool_storage(pools->feature_counts_pool);
    free_memory_pool_storage(pools->feature_umi_counts_pool);
    free_memory_pool_storage(pools->feature_sequences_pool);
    free_memory_pool_storage(pools->unmatched_barcodes_features_block_pool);
}
void free_feature_arrays(feature_arrays *features) {
    free(features->feature_names_storage);
    free(features->feature_lengths);
    free(features->feature_code_lengths);
    free(features->feature_sequences_storage);
    free(features->feature_codes_storage);
    free(features->feature_names);
    free(features->feature_sequences);
    free(features->feature_codes);
    free(features);
}

guint hash_int64(gconstpointer v) {
    guint64 k = *(const guint64*)v;
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    return k;
}
// Simple equality function for 64-bit integers
gboolean equal_int64(gconstpointer v1, gconstpointer v2) {
    return *(const long long*)v1 == *(const long long*)v2;
}

guint hash_int32(gconstpointer v) {
    guint32 k = *(const guint32*) v;
    k ^= k >> 16;
    k *= 0x85ebca6b;
    k ^= k >> 13;
    return k;
}
// Simple equality function for 64-bit integers
gboolean equal_int32(gconstpointer v1, gconstpointer v2) {
     return *((const gint*)v1) == *((const gint*)v2);
}

int find_substring_match(const char *query, feature_arrays *features) {
    for (int i = 0; i < features->number_of_features ; i++) {
       const int set_len = features->feature_lengths[i]; // Length of the current set string
        for (int j = 0; j < strlen(query); j++) {
            if (j + set_len <= strlen(query)) { // Only consider substrings as long as the set string
                char *substring = strndup(query + j, set_len); // Create a substring of length set_len
                if (strstr(features->feature_sequences[i], substring) != NULL) {
                    free(substring); // Free the dynamically allocated substring
                    return i+1; // Return the index of the set where the substring is found
                }
                free(substring); // Ensure memory is freed if not returning
            }
        }
    }
    return 0; // Return -1 if no substring match is found in any string
}


storage_block* allocate_storage_block(size_t block_size){
    storage_block *new_block = (storage_block *)malloc(sizeof(storage_block));
    if (!new_block) {
        perror("Failed to allocate new storage block");
        return NULL;
    }
    new_block->storage = (unsigned char *)malloc(block_size);
    if (!new_block->storage) {
        perror("Failed to allocate new storage block");
        return NULL;
    }
    new_block->next = NULL;
    memset(new_block->storage, 0, block_size);
    return new_block;
}
memory_pool* initialize_storage_pool(size_t block_size, size_t blocks_per_pool){
    memory_pool *pool = (memory_pool *)malloc(sizeof(memory_pool));
    if (!pool) {
        perror("Failed to allocate memory for memory pool");
        exit(EXIT_FAILURE);
    }
    pool->block_size = block_size;
    pool->blocks_per_pool = blocks_per_pool;
    pool->free_blocks = blocks_per_pool;
    pool->first_block = allocate_storage_block(block_size*blocks_per_pool);
    pool->current_block = pool->first_block;
    pool->next_free = pool->current_block->storage;
    return pool;
}

void* allocate_memory_from_pool(memory_pool *pool){
    if (pool->free_blocks == 0) {
        expand_memory_pool(pool);
    }
    pool->free_blocks--;
    void *ret = pool->next_free;
    pool->next_free += pool->block_size;
    return ret;
}
void expand_memory_pool(memory_pool *pool){
    storage_block *new_block = allocate_storage_block(pool->block_size*pool->blocks_per_pool);
    pool->current_block->next = new_block;
    pool->current_block = new_block;
    pool->next_free = pool->current_block->storage;
    pool->free_blocks += pool->blocks_per_pool;
}

void initialize_unit_sizes(){
    // Size of feature_counts (rounded up to alignment of uint16_t)
    size_t feature_counts_size = sizeof(feature_counts) + 2 * (number_of_features + 1) * sizeof(uint16_t);
    size_t feature_counts_alignment = __alignof__(uint16_t);  // Use __alignof__ instead of alignof
    dynamic_struct_sizes.feature_counts = (feature_counts_size + feature_counts_alignment - 1) & ~(feature_counts_alignment - 1);

    // Size of feature_umi_counts (rounded up to alignment of uint16_t)
    size_t feature_umi_counts_size = sizeof(feature_umi_counts) + (number_of_features + 1) * sizeof(uint16_t);
    size_t feature_umi_counts_alignment = __alignof__(uint16_t);  // Use __alignof__ instead of alignof
    dynamic_struct_sizes.feature_umi_counts = (feature_umi_counts_size + feature_umi_counts_alignment - 1) & ~(feature_umi_counts_alignment - 1);

    // Size of feature_sequences (rounded up to alignment of char)
    size_t feature_sequences_size = sizeof(feature_sequences) + maximum_feature_length + 1;
    fprintf(stderr, "Feature sequences size %ld\n", feature_sequences_size);
    size_t feature_sequences_alignment = __alignof__(char);  // Use __alignof__ instead of alignof
    dynamic_struct_sizes.feature_sequences = (feature_sequences_size + feature_sequences_alignment - 1) & ~(feature_sequences_alignment - 1);
    fprintf(stderr, "Adjusted feature sequences size %ld\n", dynamic_struct_sizes.feature_sequences);

    // Size of unmatched_barcodes_features_block (rounded up to alignment of unsigned char)
    size_t unmatched_barcodes_features_block_size = sizeof(unmatched_barcodes_features_block)
            + barcode_code_length 
            + umi_code_length
            + 2
            + (max_barcode_mismatches + 1) * barcode_code_length
            + (max_barcode_mismatches + 1);
    size_t unmatched_barcodes_features_block_alignment = __alignof__(unsigned char);  // Use __alignof__ instead of alignof
    dynamic_struct_sizes.unmatched_barcodes_features_block = (unmatched_barcodes_features_block_size + unmatched_barcodes_features_block_alignment - 1) & ~(unmatched_barcodes_features_block_alignment - 1);
}

void read_unmatched_features_block(unmatched_barcodes_features_block *entry_block, unmatched_barcodes_features *entry){
    entry->next=entry_block->next;
    entry->feature_index=entry_block->storage[0];
    entry->barcode=entry_block->storage+1;
    entry->umi=entry->barcode+barcode_code_length;
    entry->number_of_closest_barcodes=entry->umi[umi_code_length];
    entry->closest_barcodes=entry->umi+umi_code_length+1;
    entry->Qscores=entry->closest_barcodes+(barcode_code_length)*(max_barcode_mismatches+1);
}
int insert_feature_sequence(char *sequence, unsigned char feature_index, unsigned char hamming_distance, data_structures *hashes, memory_pool_collection *pools){
    gpointer *value=g_hash_table_lookup(hashes->unique_features_match, sequence);
    if (value){
        feature_sequences *entry = (feature_sequences*)value;
        entry->counts++;
        return 0;
    }
    else{
        feature_sequences *new_entry = (feature_sequences*) allocate_memory_from_pool(pools->feature_sequences_pool);
        strcpy(new_entry->sequence, sequence);
        new_entry->feature_index=feature_index;
        new_entry->hamming_distance=hamming_distance;
        new_entry->counts=1;
        g_hash_table_insert(hashes->unique_features_match, new_entry->sequence, new_entry);
        return 1;
    }   
}
int print_feature_sequences(feature_arrays *features, int *total_counts, char *directory, data_structures *hashes){
    //remember that 1 is the first feature
    memset(total_counts, 0, features->number_of_features * sizeof(int));
    char feature_sequences_filename[FILENAME_LENGTH];
    sprintf(feature_sequences_filename, "%s/feature_sequences.txt", directory);
    FILE *feature_sequencesfp = fopen(feature_sequences_filename, "w");
    if (feature_sequencesfp == NULL) {
        fprintf(stderr, "Failed to open feature sequences file\n");
        exit(EXIT_FAILURE);
    }   
    GHashTableIter iter;
    gpointer key, value;
    int i=0;
    g_hash_table_iter_init(&iter, hashes->unique_features_match);
    int count = g_hash_table_size(hashes->unique_features_match);
    gpointer *array = g_new(gpointer, count);
    if(!array){
        fprintf(stderr, "Failed to allocate memory for array\n");
        exit(EXIT_FAILURE);
    }   
    while (g_hash_table_iter_next(&iter, &key, &value)) { 
        array[i++] = value;
    }
    qsort(array, count, sizeof(gpointer), compare_feature_sequences);
    for (i = 0; i < count; i++) {
        char annotated_sequence[LINE_LENGTH];
        const int mapped_index=((feature_sequences*)array[i])->feature_index-1;
        feature_sequences *entry = (feature_sequences*)array[i];
        total_counts[mapped_index]+=entry->counts;
        strcpy(annotated_sequence, entry->sequence);
        lowerCaseDifferences(features->feature_sequences[entry->feature_index-1],annotated_sequence,strlen(annotated_sequence));
        fprintf(feature_sequencesfp, "%3d %s %2d %7d %s\n", entry->feature_index,annotated_sequence, entry->hamming_distance,entry->counts,features->feature_names[entry->feature_index-1]);
    }
    return 0;
}

int string2code_debug(char *string, int sequence_length, unsigned char *code){
    const int last_element=sequence_length/4;
    for (int i=0; i<last_element; i++){
        DEBUG_PRINT( "Converting %c %c %c %c\n", string[4*i], string[4*i+1], string[4*i+2], string[4*i+3]    );
        code[i]=seq2code[(unsigned char)string[4*i]]<<6 | seq2code[(unsigned char)string[4*i+1]]<<4 | seq2code[(unsigned char)string[4*i+2]]<<2 | seq2code[(unsigned char)string[4*i+3]];
    }
    //check if there are any remaining characters and pad with 0
    char leftover=sequence_length%4;
    if (!leftover){
        return last_element;
    }
    switch (leftover){
        case 1:
            DEBUG_PRINT( "Converting %c \n", string[4*last_element]);
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6;
            break;
        case 2:
            DEBUG_PRINT( "Converting %c %c\n", string[4*last_element], string[4*last_element+1]);
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6 | seq2code[(unsigned char)string[4*last_element+1]]<<4;
            break;
        case 3:
            DEBUG_PRINT( "Converting %c %c %c\n", string[4*last_element], string[4*last_element+1], string[4*last_element+2]);
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6 | seq2code[(unsigned char)string[4*last_element+1]]<<4 | seq2code[(unsigned char)string[4*last_element+2]]<<2;
            break;
    }
    return last_element+1;
}

unmatched_barcodes_features_block* add_unmatched_barcode_store_feature(unsigned char *barcodes, char *umi, unsigned char *qscores, unsigned char feature_index,int number_of_variants, memory_pool_collection *pools, statistics *stats){
    unmatched_barcodes_features_block *new_entry = (unmatched_barcodes_features_block*)allocate_memory_from_pool(pools->unmatched_barcodes_features_block_pool); 

    if (stats->unmatched_list.first_entry== NULL){
        stats->unmatched_list.first_entry=new_entry;
        stats->unmatched_list.last_entry=new_entry;
    }
    else{
        stats->unmatched_list.last_entry->next=new_entry;
        stats->unmatched_list.last_entry=new_entry;
    }
    unsigned char *storage = new_entry->storage; 
    storage[0]=feature_index;
    memcpy(storage+1, barcodes, barcode_code_length);
    storage+=1+barcode_code_length;
    string2code(umi, umi_length, storage);
    storage+=umi_code_length;    
    storage[0]=number_of_variants;
    storage++;
    //leave empty space for the closest barcodes    
    storage+=(barcode_code_length)*(max_barcode_mismatches+1);
    memcpy(storage, qscores, max_barcode_mismatches+1);
    new_entry->next=NULL;
    return new_entry;
}


unsigned char* readWhiteList(char *whitelist_filename,GHashTable *hash){
    FILE *whitelist_file = fopen(whitelist_filename, "r");
    if (whitelist_file == NULL) {
        fprintf(stderr, "Failed to open whitelist file %s\n", whitelist_filename);
        exit(EXIT_FAILURE);
    }
    char line[LINE_LENGTH];
    //just in case the file is corrupted - calculate the sizes by counting lines
    //find the number of lines in the file by counting the number of newlines
    char ch=0,lastChar=0;
    size_t line_count=0;
    //from the first line calculate the barcode length and check if the sequence is valid
    int sequence_count=0;
    while((ch = fgetc(whitelist_file)) != EOF) {
        if(ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'N'){
            sequence_count++;
        }
        else if (ch == '\n') {
            if(sequence_count){
                if (barcode_length==0){
                    barcode_length=sequence_count;
                    barcode_code_length=(barcode_length+3)/4;
                }
                else if (sequence_count != barcode_length){
                    fprintf(stderr, "Error: Invalid barcode length %d\n", sequence_count);
                    exit(EXIT_FAILURE);
                }
                sequence_count=0;
                line_count++;
            }
        }
        else{
            fprintf(stderr, "Error: Invalid character %c in whitelist file\n", ch);
            exit(EXIT_FAILURE);
        }
    }
    // If the last character isn't a newline, increment the count
    if (lastChar != '\n') {
        line_count++;
    }
    
    whitelist = malloc(line_count * barcode_code_length);
    if (whitelist == NULL) {
        fprintf(stderr, "Failed to allocate memory for whitelist storage\n");
        exit(EXIT_FAILURE);
    }
    memset(whitelist, 0, line_count * barcode_code_length);
    //reset the file pointer to the beginning of the file
    fseek(whitelist_file, 0, SEEK_SET);
    size_t nBarcodes=0;
    //set line to zero
    memset(line, 0, LINE_LENGTH);
    while ( fgets(line, LINE_LENGTH, whitelist_file) != NULL) {
        if (!check_sequence(line, barcode_length)){
            fprintf(stderr, "Error: Invalid barcode sequence of expected length %d %s\n",barcode_length,line);
            exit(EXIT_FAILURE);
        }
        //up to 16 characters will be stored in uint32_t
        memset(line+barcode_length, 0, LINE_LENGTH-barcode_length);
        int j=0;
        int i=0;
        const size_t offset=nBarcodes*barcode_code_length;
        memset(whitelist+offset, 0, barcode_code_length);
        while(j<barcode_length){
            unsigned char char_value=seq2code[(unsigned char)line[j]]<<6 | seq2code[(unsigned char)line[j+1]]<<4 | seq2code[(unsigned char)line[j+2]]<<2 | seq2code[(unsigned char)line[j+3]];
            whitelist[offset+i]=char_value;
            i++;
            j+=4;
        }
        if(!g_hash_table_insert(hash,  (uint32_t*)(whitelist+offset) , whitelist+offset)){
            fprintf(stderr, "Error: Failed to insert barcode %s into the hash table %ld\n", line, nBarcodes);
            exit(EXIT_FAILURE); 
        }
        nBarcodes++;
    }
    fprintf(stderr, "Read %ld barcodes\n", nBarcodes);
    fclose(whitelist_file);
    return whitelist;
}

feature_arrays* readTagsFile(const char* filename) {
    int seq_size=0;
    int name_size=0;
    int code_size=0;
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open tags file");
        exit(EXIT_FAILURE);
    }
    char line[LINE_LENGTH];
    int count = 0;
    //skip the header
    //count the lines and check that the sequences are valid
    int maxFeatureLength=0;
    if (!fgets(line, LINE_LENGTH, file)) {
        perror("Failed to read tags file");
        exit(EXIT_FAILURE);
    }
    while (fgets(line, LINE_LENGTH, file) != NULL) {
        // Split the line by spaces and read the 3rd and 6th columns
        char *token;
        int colIndex = 0;
        char *tmpSeq = NULL;

        token = strtok(line, ","); // Assuming comma-separated values as per the problem statement
        while (token != NULL) {
            colIndex++;
            if (colIndex == 2) { // Second column for names
                name_size+=strlen(token)+1;
            } 
            else if (colIndex == 5) { // Third column for sequences
                tmpSeq = token;
                // Remove possible newline character from the sequence
                if (!check_sequence(tmpSeq, strlen(tmpSeq))){
                    //fprintf(stderr, "Error: Invalid sequence %s\n", tmpSeq);
                    exit(EXIT_FAILURE);
                }
                const int string_length=strlen(tmpSeq);
                seq_size+=strlen(tmpSeq)+1; 
                code_size+=string_length/4;
                if (string_length%4){
                    code_size++;
                } 
                if (strlen(tmpSeq) > maxFeatureLength){
                    maxFeatureLength=strlen(tmpSeq);
                }
                break;
            }
            token = strtok(NULL, ",");
        }
        count++;
    }
    fprintf(stderr, "Read %d tags with max length %d\n", count, maxFeatureLength);
    feature_arrays *myfeatures = malloc(sizeof(feature_arrays));
    if (myfeatures == NULL) {
        fprintf(stderr, "Failed to allocate memory for feature arrays\n");
        exit(EXIT_FAILURE);
    }
    memset(myfeatures, 0, sizeof(feature_arrays));
    myfeatures->max_length=maxFeatureLength;
    myfeatures->feature_names_storage=malloc(name_size);
    myfeatures->feature_sequences_storage=malloc(seq_size);
    myfeatures->feature_codes_storage=malloc(code_size);
    myfeatures->feature_names=malloc(count*sizeof(char*));
    myfeatures->feature_lengths=malloc(count*sizeof(unsigned int));
    myfeatures->feature_code_lengths=malloc(count*sizeof(unsigned char));
    myfeatures->feature_codes=malloc(count*sizeof(unsigned char*));
    myfeatures->feature_sequences=malloc(count*sizeof(char*));
    myfeatures->number_of_features=count;

    //check if any of the mallocs failed by checking for 0 pointers
    if (myfeatures->feature_names_storage == NULL || myfeatures->feature_sequences_storage == NULL || myfeatures->feature_codes_storage == NULL || myfeatures->feature_names == NULL || myfeatures->feature_lengths == NULL || myfeatures->feature_code_lengths == NULL || myfeatures->feature_codes == NULL){
        fprintf(stderr, "Failed to allocate memory for feature arrays\n");
        exit(EXIT_FAILURE);
    }
    memset(myfeatures->feature_names_storage, 0, name_size);
    memset(myfeatures->feature_sequences_storage, 0, seq_size);
    memset(myfeatures->feature_codes_storage, 0, code_size);
    myfeatures->feature_names[0]=myfeatures->feature_names_storage;
    myfeatures->feature_sequences[0]=myfeatures->feature_sequences_storage;
    myfeatures->feature_codes[0]=myfeatures->feature_codes_storage;

    //rewind the file and read the sequences
    fseek(file, 0, SEEK_SET);
    if (!fgets(line, LINE_LENGTH, file)) {
        perror("Failed to read tags file");
        exit(EXIT_FAILURE);
    }
    count=0;
    while (fgets(line, LINE_LENGTH, file) != NULL) {
        // Split the line by spaces and read the 3rd and 6th columns
        char *token;
        int colIndex = 0;
        char *tmpName = NULL;
        char *tmpSeq = NULL;

        token = strtok(line, ","); // Assuming comma-separated values as per the problem statement
        while (token != NULL) {
            colIndex++;
            if (colIndex == 2) { // Second column for names
                tmpName = token;
                strcpy(myfeatures->feature_names[count],tmpName);
                if(count + 1 < myfeatures->number_of_features){
                    myfeatures->feature_names[count+1]=myfeatures->feature_names[count]+strlen(tmpName)+1;
                }
            } else if (colIndex == 5) { // Third column for sequences
                tmpSeq = token;
                strcpy(myfeatures->feature_sequences[count],tmpSeq);
                myfeatures->feature_lengths[count]=strlen(tmpSeq);
                myfeatures->feature_code_lengths[count]=string2code(tmpSeq, strlen(tmpSeq),myfeatures->feature_codes[count]);
                if(count + 1 < myfeatures->number_of_features){
                    myfeatures->feature_sequences[count+1]=myfeatures->feature_sequences[count]+strlen(tmpSeq)+1;
                    myfeatures->feature_codes[count+1]=myfeatures->feature_codes[count]+myfeatures->feature_code_lengths[count];
                }
                break;
            }
            token = strtok(NULL, ",");
        }
        count++;
    }
    fclose(file);
    return myfeatures;
}
char check_sequence(char *sequence, int sequence_length){
    for (int i=0; i<sequence_length; i++){
        if (sequence[i] == 'A' || sequence[i] == 'C' || sequence[i] == 'G' || sequence[i] == 'T'){
            continue;
        }
        return 0;
    }
    return 1;
}
int string2code(char *string, int sequence_length, unsigned char *code){
    const int last_element=sequence_length/4;
    for (int i=0; i<last_element; i++){
        code[i]=seq2code[(unsigned char)string[4*i]]<<6 | seq2code[(unsigned char)string[4*i+1]]<<4 | seq2code[(unsigned char)string[4*i+2]]<<2 | seq2code[(unsigned char)string[4*i+3]];
    }
    //check if there are any remaining characters and pad with 0
    char leftover=sequence_length%4;
    if (!leftover){
        return last_element;
    }
    switch (leftover){
        case 1:
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6;
            break;
        case 2:
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6 | seq2code[(unsigned char)string[4*last_element+1]]<<4;
            break;
        case 3:
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6 | seq2code[(unsigned char)string[4*last_element+1]]<<4 | seq2code[(unsigned char)string[4*last_element+2]]<<2;
            break;
    }
    return last_element+1;
}
void string2all_codes(char *string, unsigned char codes[][LINE_LENGTH/2+1], int *lengths){
    //4 codes are returned for the string to capture all the possible frames
    char offset[3][LINE_LENGTH];
    const int seqlength = strlen(string);

    strcpy(offset[0], string+1);
    strcpy(offset[1], string+2);
    strcpy(offset[2], string+3);
    lengths[0]=string2code(string, seqlength, codes[0]);

    for (int i=0; i<3; i++){
        lengths[i+1]=string2code(offset[i], seqlength-1-i, codes[i+1]);
    }   
}

int check_barcode(char *barcode, unsigned char *code){
    for (int i=0; i<barcode_code_length; i++){
        string2code(barcode, barcode_length, code);
    }
    if (g_hash_table_lookup(whitelist_hash,code) == NULL){
        return 1;
    }
    return 0;
}
int find_matching_codes(unsigned char *sequence_code, unsigned char *feature_code, size_t feature_code_length, size_t feature_length){
    unsigned char a1[feature_code_length];
    unsigned char *b=feature_code;
    memcpy(a1, sequence_code, feature_code_length);
    unsigned char right_overhang=feature_length % 4;

    if (right_overhang){
        a1[feature_code_length-1] = a1[feature_code_length-1] & (0xff << (8-2*right_overhang));        
    }
    int length=feature_code_length;
    unsigned char *a=a1;
    if(length >= 8){
        uint64_t *a64=(uint64_t*) a;
        uint64_t *b64=(uint64_t*) b;
        for (int i=0; i<length/8; i++){
            if (a64[i] != b64[i]){
                return feature_code_length+1;
            }
        }
        if (length % 8 == 0){
            return 0;
        }
        a+=8*(length/8);
        b+=8*(length/8);
        length=length%8;
    }
    if(length >= 4){
        const uint32_t *a32=(uint32_t*) a;
        const uint32_t *b32=(uint32_t*) b;
        for (int i=0; i<length/4; i++){
            if (a32[i] != b32[i]){
                return feature_code_length+1;
            }
        }
        if (length % 4 == 0){
            return 0;
        }
        a+=4*(length/4);
        b+=4*(length/4);
        length=length%4;    
    }
    if (length >= 2){
        const uint16_t *a16=(uint16_t*) a;
        const uint16_t *b16=(uint16_t*) b;
        for (int i=0;i<length/2; i++){
            if (a16[i] != b16[i]){
                return feature_code_length+1;  
            }
        }
        if (length % 2 == 0){
            return 0;
        }
        a+=2*(length/2);
        b+=2*(length/2);
        length=length%2;
    }
    for (int i=0;i<length; i++){
        if (a[i] != b[i]){
            return length+1;
        }
    }
    return 0;
}
int fuzzy_matching_codes(unsigned char *sequence_code, unsigned char *feature_code, size_t feature_code_length, size_t feature_length, int maxHammingDistance){
    int mismatches=0;
    unsigned char a1[feature_code_length];
    unsigned char *b=feature_code;
    memcpy(a1, sequence_code, feature_code_length);
    unsigned char right_overhang=feature_length % 4;
    if (right_overhang){
        a1[feature_code_length-1] = a1[feature_code_length-1] & (0xff << (8-2*right_overhang));        
    }
    unsigned char *a=a1;
    int length=feature_code_length;
    if(length >= 8){
        const uint64_t *a64=(uint64_t*) a;
        const uint64_t *b64=(uint64_t*) b;
        for (int i=0; i<length/8; i++){
            if (a64[i] != b64[i]){
                if(mismatches==maxHammingDistance){
                    return feature_code_length+1;
                }
                uint64_t diff=a64[i]^b64[i];
                unsigned char *diff8=(unsigned char*)&diff;
                mismatches+=diff2Hamming[diff8[0]]+diff2Hamming[diff8[1]]+diff2Hamming[diff8[2]]+diff2Hamming[diff8[3]] + diff2Hamming[diff8[4]]+diff2Hamming[diff8[5]]+diff2Hamming[diff8[6]]+diff2Hamming[diff8[7]];
                if (mismatches > maxHammingDistance){
                    return feature_code_length+1;
                }
            }
        }
        if (length % 8 == 0){
            return mismatches;
        }
        a+=8*(length/8);
        b+=8*(length/8);
        length=length%8;
    }
    if (length >= 4 ){
        const int last=8*(length/8);
        const uint32_t *a32=(uint32_t*) a;
        const uint32_t *b32=(uint32_t*) b;
        for (int i=last/4; i<length/4; i++){
            if (a32[i] != b32[i]){
                if(mismatches==maxHammingDistance){
                    return feature_code_length+1;
                }
                uint32_t diff=a32[i]^b32[i];
                unsigned char *diff8=(unsigned char*)&diff;
                mismatches+=diff2Hamming[diff8[0]]+diff2Hamming[diff8[1]]+diff2Hamming[diff8[2]]+diff2Hamming[diff8[3]];
                if (mismatches > maxHammingDistance){
                    return feature_code_length+1;
                }
            }
        }
        if (length % 4 == 0){
            return mismatches;
        }
        a+=4*(length/4);
        b+=4*(length/4);
        length=length%4;
        
    }
    if (length >= 2){
        const uint16_t *a16=(uint16_t*) a;
        const uint16_t *b16=(uint16_t*) b;
        const int last=4*(length/4);
        for (int i=last/2;i<length/2; i++){
            if (a16[i] != b16[i]){
                if(mismatches==maxHammingDistance){
                    return feature_code_length+1;
                }
                uint16_t diff=a16[i]^b16[i];
                unsigned char *diff8=(unsigned char*)&diff;
                mismatches+=diff2Hamming[diff8[0]]+diff2Hamming[diff8[1]];  
                if (mismatches > maxHammingDistance){
                    return feature_code_length+1;
                }   
            }
        }
        if (length % 2 == 0){
            return mismatches;
        }
        a+=2*(length/2);
        b+=2*(length/2);
        length=length%2;
    }
    const int last=2*(length/2);
    for (int i=last;i<length; i++){
        if (a[i] != b[i]){
            if(mismatches==maxHammingDistance){
                return feature_code_length+1;
            }
            mismatches+=diff2Hamming[a[i]^b[i]];
            if (mismatches > maxHammingDistance){
                return feature_code_length+1;
            }
        }
    }   
    return mismatches;
}

int find_matches_in_sub_arrays(unsigned char *sequence_code, unsigned char *feature_code, size_t sequence_code_length, size_t feature_code_length, size_t feature_length, int maxHammingDistance, int *best_offset){
    int offset=0;
    int length= (sequence_code_length > feature_code_length) ? sequence_code_length : feature_code_length;
    int minHammingDistance=length+1;
    *best_offset=0;
    while (sequence_code_length - offset >= feature_code_length){
        if (!maxHammingDistance){
            if (!find_matching_codes(sequence_code+offset, feature_code, feature_code_length, feature_length)){
                *best_offset=offset;
                return 0;
            }
        }
        else{
            int hammingDistance=fuzzy_matching_codes(sequence_code+offset, feature_code, feature_code_length, feature_length, maxHammingDistance);
            if (hammingDistance < minHammingDistance){
                minHammingDistance=hammingDistance;
                *best_offset=offset;
            }
            //if we have found a perfect match return immediately
            if  (!minHammingDistance){
                *best_offset=offset;
                return 0;
            }   
        }
        offset++;
    }
    return minHammingDistance;
}

int checkSequenceAndCorrectForN(char *line, char *corrected_lines[], char *buffer,int sequence_length, int maxN){
    int nCount=0;
    int indices[maxN];
    corrected_lines[0]=line;
    for (int i=0; i<sequence_length; i++){
        if (line[i] == 'N'){
            indices[nCount++]=i;
            if (nCount >=  maxN) return 0;
            
        }
        else if (line[i] != 'A' && line[i] != 'C' && line[i] != 'G' && line[i] != 'T'){
            return 0;
        }
    }
    if (!nCount){
        return 1;
    }
    char line_copy[LINE_LENGTH];
    strcpy(line_copy, line);
    return generate_sequences(line_copy, sequence_length, indices, &buffer, 0, nCount, corrected_lines, 0);
}char* printToBuffer(char *string, int sequence_length, char *buffer){
    for (int i=0; i<sequence_length; i++){
        buffer[i]=string[i];
    }
    buffer[sequence_length]='\0';
    return buffer + sequence_length+1;
}

int generate_sequences(char* string, int string_length, int *indices, char **output, int index, int number_of_indices, char *output_indices[], 
int number_of_outputs) {
    static char bases[] = "ACGT";
    //base case when index is equal to number of indices
    if (index == number_of_indices) {
        printToBuffer(string, string_length, *output);
        output_indices[number_of_outputs++]=*output;
        *output += string_length + 1;
        return number_of_outputs;
    }
    //when the index is less than the number of indices
    for (int i = 0; i < 4; i++) {
        string[indices[index]] = bases[i];  
        number_of_outputs=generate_sequences(string, string_length, indices, output, index+1, number_of_indices, output_indices, number_of_outputs);
    }
    return number_of_outputs;
}
int simple_search(feature_arrays *features, char *line){  
    for (int i=0; i<features->number_of_features; i++){
        char *query=line;
        char *feature=features->feature_sequences[i];
        //fprintf(stderr, "Comparing %s %s\n", query, feature);
        while (*query == *feature){
            query++;
            feature++;
        }
        if (!*feature){
            return i+1;
        }
    }
    return 0;
}
int simple_hamming_search(feature_arrays *features, char *line, int maxHammingDistance, int *hamming_distance){  
    int ambiguous=0;
    int bestFeature=0;
    int bestHammingDistance=maxHammingDistance+1;
    for (int i=0; i<features->number_of_features; i++){
        char *query=line;
        char *feature=features->feature_sequences[i];
        int hammingDistance=0;
        
        //fprintf(stderr, "Comparing %s %s\n", query, feature);
        while (*query && *feature){
            if (*query != *feature){
                hammingDistance++;
                if (hammingDistance > maxHammingDistance){
                    break;
                }
            }
            query++;
            feature++;
        }
        if (!*feature){
           if (hammingDistance == bestHammingDistance){
               ambiguous=1;
           }
           else if (hammingDistance < bestHammingDistance){
               bestFeature=i+1;
               bestHammingDistance=hammingDistance;
               ambiguous=0;
           }
        }

    }
    *hamming_distance=bestHammingDistance;
    if (bestHammingDistance <= maxHammingDistance && !ambiguous){
        return bestFeature;
    }
    return 0;
}
int find_feature_match_single(feature_arrays *features, char *lineR2, int maxHammingDistance,int *bestScore, char **matching_sequence){
    // convert lineR2 to 4 codes
    // do a quick check to see if there is a perfect match the constant feature
    //int best_feature=simple_search_code(features, lineR2, constant_offset);
    unsigned char codes[4][LINE_LENGTH/2+1];
    int code_lengths[4];
    string2all_codes(lineR2, codes, code_lengths);
    int best_feature=0;
    int bestHammingDistance=maxHammingDistance+1; 
    int best_sequence_offset=0;
    int ambiguous=0;
    for (int i=0; i<4; i++){
        for (int j=0; j<features->number_of_features; j++){
            int code_offset=0;
            int hammingDistance=find_matches_in_sub_arrays(codes[i], features->feature_codes[j], code_lengths[i], features->feature_code_lengths[j],features->feature_lengths[j], maxHammingDistance, &code_offset);
            if (!hammingDistance){
                *bestScore=0;
                best_feature=j+1;
                *matching_sequence=lineR2+i+code_offset*4;
                return best_feature;
            }
            if (hammingDistance < bestHammingDistance){
                best_feature=j+1;
                bestHammingDistance=hammingDistance;
                best_sequence_offset=i+4*code_offset;
                ambiguous=0;
            }
            else if (hammingDistance == bestHammingDistance){
                ambiguous=1;
            }
        }
    }
    if (bestHammingDistance <= maxHammingDistance){
        *bestScore=bestHammingDistance;
        *matching_sequence=lineR2+best_sequence_offset;
        if (!ambiguous){
            return best_feature;
        }
        return 0;
    }
    return 0;
}
int find_feature_match_parallel(feature_arrays *features, char *lineR2, int maxHammingDistance, int nThreads, int *bestScore, char **matching_sequence){
    // convert lineR2 to 4 codes
    // do a quick check to see if there is a perfect match the constant feature
    //int best_feature=simple_search_code(features, lineR2, constant_offset);
    nThreads=(nThreads > 4) ? 4 : nThreads;
    nThreads=(nThreads < 1) ? 1 : nThreads;
    if (nThreads==1){
        return find_feature_match_single(features, lineR2, maxHammingDistance, bestScore, matching_sequence);
    }
    unsigned char codes[4][LINE_LENGTH/2+1];
    int code_lengths[4];
    string2all_codes(lineR2, codes, code_lengths);
    int best_feature=0;
    int bestFeatureDistance=maxHammingDistance+1;   
    int bestHammingDistances[4]={maxHammingDistance+1,maxHammingDistance+1,maxHammingDistance+1,maxHammingDistance+1};
    int best_code_offsets[4]={0,0,0,0};
    int ambiguous[4]={0,0,0,0};
    int best_match[4]={0,0,0,0};

    int exact_match_found=0;
    #pragma omp parallel for num_threads(nThreads)
    for (int i=0; i<4; i++){
        //get the thread number
        const int thread_num=omp_get_thread_num();
        for (int j=0; j<features->number_of_features && !exact_match_found; j++){
            int code_offset=0;
            int hammingDistance=find_matches_in_sub_arrays(codes[i], features->feature_codes[j], code_lengths[i], features->feature_code_lengths[j],features->feature_lengths[j], maxHammingDistance, &code_offset);
            //fprintf(stderr, "Thread %d feature %d hamming distance %d\n", thread_num, j, hammingDistance);
            if (!hammingDistance){
                *bestScore=0;
                exact_match_found=1;
                best_feature=j+1;
                *matching_sequence=lineR2+i+code_offset*4;
                break;
            }
            if (hammingDistance < bestHammingDistances[thread_num]){
                best_match[thread_num]=j+1;
                bestHammingDistances[thread_num]=hammingDistance;
                best_code_offsets[thread_num]=code_offset;
                ambiguous[thread_num]=0;
            }
            else if (hammingDistance == bestHammingDistances[thread_num]){
                ambiguous[thread_num]=1;
            }
        }
    }
    if (exact_match_found){
        return best_feature;
    }
    //find the best hamming distance - if there is a perfect match return immediately
    char multiAmbiguous=0;
    int best_code_offset=0;
    int best_i=0;
    for (int i=0; i<nThreads; i++){
        if (!bestHammingDistances[i]){
            *bestScore=0;
            *matching_sequence=lineR2+i+best_code_offsets[i]*4;
            return best_match[i];   
        }
        if (bestHammingDistances[i] < bestFeatureDistance){
            bestFeatureDistance=bestHammingDistances[i];
            if(ambiguous[i]){
                multiAmbiguous=1;
            }
            else{
                best_feature=best_match[i];
                best_code_offset=best_code_offsets[i];
                multiAmbiguous=0; //reset ambiguous flag
                best_i=i;
            }
        }
        else if (bestHammingDistances[i] == bestFeatureDistance){
            multiAmbiguous=1;
        }
    }
    if (bestFeatureDistance <= maxHammingDistance){
        *bestScore=bestFeatureDistance;
        *matching_sequence=lineR2+best_i+best_code_offset*4;
        if (!multiAmbiguous){
            return best_feature;
        }
        return 0;
    }
    return 0;
}
void update_feature_counts(char *barcodeString, char *umi, uint32_t feature_index, data_structures *hashes,  memory_pool_collection *pools){ 
    unsigned char code[barcode_code_length];
    string2code(barcodeString, barcode_length, code);
    update_feature_counts_from_code(code, umi, feature_index, hashes, pools);
}
void update_feature_counts_from_code(unsigned char *code, char *umi,  uint32_t feature_index, data_structures *hashes, memory_pool_collection *pools){
    update_umi_counts(code, umi, feature_index, hashes, pools);
    feature_counts *s=g_hash_table_lookup(hashes->filtered_hash, code);
    if (s == NULL) {
        feature_counts *entry= (feature_counts*) allocate_memory_from_pool(pools->feature_counts_pool);
        memcpy(entry->sequence_code,code, barcode_code_length );
        entry->counts[feature_index]=1;
        g_hash_table_insert(hashes->filtered_hash, entry->sequence_code, entry);
        return;
    }    
    s->counts[feature_index]++;
}
void update_umi_counts(unsigned char *code, char *umi,  uint32_t feature_index,data_structures *hashes, memory_pool_collection *pools){
    unsigned char code8[8];
    memset(code8,0,8);
    memcpy(code8,code,barcode_code_length);
    string2code(umi, umi_length, code8+barcode_code_length);
    feature_umi_counts *s=g_hash_table_lookup(hashes->sequence_umi_hash, code8);
    if(s == NULL){
        feature_umi_counts *entry= (feature_umi_counts*) allocate_memory_from_pool(pools->feature_umi_counts_pool);
        memset(entry,0,sizeof(feature_umi_counts));
        memcpy(entry->sequence_umi_code,code8, 8);
        entry->counts[feature_index]=1;
        g_hash_table_insert(hashes->sequence_umi_hash, entry->sequence_umi_code, entry);
        return;
    }
    s->counts[feature_index]++;
}
char check_neighbor(uint64_t code64,uint32_t *counts, data_structures *hashes){
    feature_umi_counts *result=g_hash_table_lookup(hashes->sequence_umi_hash, &code64);
    if (result && !result->counts[0]){
        for (int i=0; i<number_of_features+1; i++){
            counts[i]+=result->counts[i];
        }
        //mark as visited
        result->counts[0]=1;
        return 1;
    }
    return 0;
}
int find_neighbors(uint64_t key64, uint64_t *neighbors, uint32_t *counts, data_structures *hashes){
    int neighbor_count=0;
    uint64_t code64=key64;
    unsigned char *code8= (unsigned char*) &code64;
    unsigned char *code=code8+barcode_code_length;
    for (int i=0; i<umi_code_length; i++){       
        const unsigned char mask=code[i];
        for (unsigned char j=1; j<4; j++){
            code[i]=mask ^ j;
            if (check_neighbor(code64,counts,hashes)) neighbors[neighbor_count++]=code64;
            code[i]=mask ^ (j << 2);
            if (check_neighbor(code64,counts,hashes)) neighbors[neighbor_count++]=code64;
            code[i]=mask ^ (j << 4);
            if (check_neighbor(code64,counts,hashes)) neighbors[neighbor_count++]=code64;
            code[i]=mask;
        }
    }
    return neighbor_count;
}
void find_deduped_counts(data_structures *hashes ){
    GHashTableIter iter;
    gpointer lookup_key, result;
    g_hash_table_iter_init(&iter, hashes->sequence_umi_hash);
    uint32_t counts[number_of_features+1];
    while (g_hash_table_iter_next(&iter, &lookup_key, &result)) {
        if (result){
            feature_umi_counts *umi_counts = (feature_umi_counts*) result;
            if (umi_counts->counts[0]){
                continue;
            }
            memset(counts,0,sizeof(counts));
            find_connected_component(lookup_key, counts, hashes);
            feature_counts *feature_counts = g_hash_table_lookup(hashes->filtered_hash, umi_counts->sequence_umi_code);
            add_deduped_count(feature_counts, counts, stringency, min_counts);
        }
    }
}

void find_connected_component(gpointer start_key, uint32_t *counts, data_structures *hashes){
    uint64_t neighbors[umi_length*3];
    clear_queue(hashes->neighbors_queue);
    enqueue(hashes->neighbors_queue,*(uint64_t*)start_key);
    while (!is_empty(hashes->neighbors_queue)) {
        uint64_t code = dequeue(hashes->neighbors_queue);
        check_neighbor(code,counts,hashes);
        int neighbor_count=find_neighbors(code,neighbors, counts, hashes);
        //if (neighbor_count) fprintf(stderr, "Found %d neighbors\n",neighbor_count);
        for (int i = 0; i < neighbor_count; i++) {
            enqueue(hashes->neighbors_queue,neighbors[i]);
        }
    }
}
void  add_deduped_count(feature_counts *s, uint32_t *counts, uint16_t stringency, uint32_t min_counts){
    uint32_t total_counts=0;
    uint16_t *const deduped_counts=&(s->counts[number_of_features+1]);
    if(!stringency){
        //all the counts are summed regardless

        for (int i=1; i<number_of_features+1; i++){
            if (counts[i] > min_counts )deduped_counts[i]++;
        }
        return;
    }
    if (stringency > 999){
        //only add if there are no duplicates
        uint32_t feature_index=0;
        total_counts=0;
        for (int i=1; i<number_of_features+1; i++){
            if (counts[i]){
                total_counts+=counts[i];    
                if (feature_index){
                    return;
                }
                feature_index=i;
            }
        }
        if (feature_index && total_counts > min_counts){
            deduped_counts[feature_index]++;
        }
        return;
    }
    //stringency is between 1 and 999
    uint32_t feature_index=0;
    
    uint32_t max_counts=0;
    unsigned char unique=0;
    for (int i=1; i<number_of_features+1; i++){
        if (counts[i]){
            total_counts+=counts[i];
            if (counts[i] > max_counts){
                max_counts=counts[i];
                feature_index=i;
                unique=1;
            }
            else if (counts[i] == max_counts){
                unique=0;
            }
        }
    }
    if (unique && (double) max_counts > total_counts*stringency/(double)1000 && total_counts > min_counts){
        deduped_counts[feature_index]++;
    }
    return;  
}

void code2string(unsigned char *code, char *string, int length){
    for (int i=0; i<length; i++){
        string[4*i]=code2seq[code[i]][0];
        string[4*i+1]=code2seq[code[i]][1];
        string[4*i+2]=code2seq[code[i]][2];
        string[4*i+3]=code2seq[code[i]][3];
    }
    string[4*length]='\0';
}

void printFeatureCounts(feature_arrays *features, int *deduped_counts, int *barcoded_counts,int **coexpression_counts, int **coexpression_histograms, char *directory, data_structures *hashes, statistics *stats){
    int total_barcodes=0;
    int total_feature_counts=0;
    int total_deduped_counts=0;
    memset(deduped_counts, 0, features->number_of_features*sizeof(int));
    memset(barcoded_counts, 0, features->number_of_features*sizeof(int));
    //print out the number of keys in filtered hash
    char barcodes_file[FILENAME_LENGTH];
    char stats_file[FILENAME_LENGTH];
    mkdir_p(directory);
    sprintf(barcodes_file, "%s/barcodes.txt", directory);
    sprintf(stats_file, "%s/stats.txt", directory);
    FILE *barcodesfp = fopen(barcodes_file, "w");
    if (barcodesfp == NULL) {
        perror("Failed to open barcodes file");
        exit(EXIT_FAILURE);
    }
    char matrix_file[LINE_LENGTH];
    sprintf(matrix_file, "%s/features_matrix.mtx", directory);
    FILE *matrixfp = fopen(matrix_file, "w");
    if (matrixfp == NULL) {
        perror("Failed to open matrix file");
        exit(EXIT_FAILURE);
    }
    find_deduped_counts(hashes);
    //print the barcodes and counts
    //iterate through the hash table
    GHashTableIter iter;
    gpointer key, value;
    g_hash_table_iter_init(&iter, hashes->filtered_hash);
    fprintf(matrixfp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(matrixfp, "%%metadata_json: {\"software_version\": \"assignBarcodes-0.1\", \"format_version\": 1}");
    
    size_t number_of_features_seen=0;
    size_t number_of_barcode_counts=0;
    while (g_hash_table_iter_next(&iter, &key, &value)) {
        feature_counts *entry = (feature_counts*)value;
        uint16_t *const entry_deduped_counts=&(entry->counts[number_of_features+1]);
        int number_of_local_features_seen=0;
        int coCounts[features->number_of_features+1];
        int coExpressorsIndices[features->number_of_features];
        memset(coCounts, 0, (features->number_of_features+1)*sizeof(int));
        int number_of_coexpressors=0;
        for (int i=1; i<= features->number_of_features; i++){
            if (entry_deduped_counts[i] > 0){
                coCounts[number_of_coexpressors]=entry_deduped_counts[i];
                coExpressorsIndices[number_of_coexpressors++]=i;
                number_of_features_seen++;
                number_of_local_features_seen++;
            }
        }
        if (number_of_local_features_seen){
            number_of_barcode_counts++;
        }
        if (number_of_coexpressors >1){
            for (int i=0; i<number_of_coexpressors; i++){
                for (int j=0; j<number_of_coexpressors; j++){
                    coexpression_counts[coExpressorsIndices[i]][coExpressorsIndices[j]]+=coCounts[i];
                }
                coexpression_histograms[coExpressorsIndices[i]][number_of_coexpressors]++;
            }
        }
        else if (number_of_coexpressors){
            coexpression_counts[coExpressorsIndices[0]][0]+=coCounts[0];
            coexpression_counts[0][coExpressorsIndices[0]]+=coCounts[0];
            coexpression_histograms[coExpressorsIndices[0]][1]++;
        }
    }
    fprintf(matrixfp, "\n%d %ld %ld\n", features->number_of_features, number_of_barcode_counts, number_of_features_seen);
    int line_no=1; 
    g_hash_table_iter_init(&iter, hashes->filtered_hash);
    while (g_hash_table_iter_next(&iter, &key, &value)) {
        feature_counts *entry = (feature_counts*)value;
        uint16_t *const entry_deduped_counts=&(entry->counts[number_of_features+1]);
        //print the barcode and counts to the matrix file
        int number_of_local_features_seen=0;
        stats->total_unmatched_features+=entry->counts[0];
        //remember that we start at feature 1
        int feature_counts[features->number_of_features+1];
        memset(feature_counts, 0, (features->number_of_features+1)*sizeof(int));
        int feature_deduped_counts[features->number_of_features+1];
        memset(feature_deduped_counts, 0, (features->number_of_features+1)*sizeof(int));
        for (int i=1; i<= features->number_of_features; i++){
            if (entry->counts[i] > 0){
                feature_counts[i]+=entry->counts[i];
                feature_deduped_counts[i]+=entry_deduped_counts[i];
            }
        }
        for (int i=1; i<= features->number_of_features; i++){
            if (feature_counts[i] > 0){
                if (feature_deduped_counts[i] > 0){
                    fprintf(matrixfp, "%d %d %d\n",i,line_no, feature_deduped_counts[i]);
                    number_of_local_features_seen++;
                }
                total_feature_counts+=feature_counts[i];
                total_deduped_counts+=feature_deduped_counts[i];
                deduped_counts[i-1]+=feature_deduped_counts[i];
                barcoded_counts[i-1]+=feature_counts[i];
            }
        }
        if (number_of_local_features_seen){
            total_barcodes++;
            line_no++;
            //print the barcode out to barcodes.txt
            char barcode[barcode_length+1];
            code2string(entry->sequence_code, barcode, barcode_code_length );
            fprintf(barcodesfp, "%s\n", barcode);
        }
    }
        
    fclose(barcodesfp);
    fclose(matrixfp);
    FILE *statsfp = fopen(stats_file, "w");
    fprintf (stderr, "Total feature counts %d\n", total_feature_counts);
    fprintf (stderr, "Total deduped feature counts %d\n", total_deduped_counts);
    fprintf (stderr, "Total unique barcode UMIs %d\n", g_hash_table_size(hashes->sequence_umi_hash));
    fprintf (stderr, "Total whitelisted barcodes %d\n", g_hash_table_size(hashes->filtered_hash));
    fprintf (stderr, "Percentage features assigned to barcode %.2f\n", 100.0*(total_feature_counts-stats->total_unmatched_features)/(double) total_feature_counts);
    fprintf (stderr, "Percentage barcodes assigned feature %.2f\n", 100.0*total_barcodes/(double) g_hash_table_size(hashes->filtered_hash));        
    fprintf (statsfp, "Total feature counts %d\n", total_feature_counts);
    fprintf (statsfp, "Total deduped feature counts %d\n", total_deduped_counts);
    fprintf (statsfp, "Total unique barcode UMIs %d\n", g_hash_table_size(hashes->sequence_umi_hash));
    fprintf (statsfp, "Total whitelisted barcodes %d\n", g_hash_table_size(hashes->filtered_hash));
    fprintf (statsfp, "Percentage features assigned to barcode %.2f\n", 100.0*(total_feature_counts-stats->total_unmatched_features)/(double) total_feature_counts);
    fprintf (statsfp, "Percentage barcodes assigned feature %.2f\n", 100.0*total_barcodes/(double) g_hash_table_size(hashes->filtered_hash));
    fclose(statsfp);

}
int find_closest_barcodes(unsigned char* code,unsigned char *corrected_codes, unsigned char *indices){
    int number_of_variants=0;
    for (int i=0; i<barcode_code_length; i++){
        unsigned char corrected_bases[3];
        //check whitelist for variant match
        int nmatches=find_variant_match(code, i, corrected_bases);
        if(nmatches){
            if (nmatches+number_of_variants > max_barcode_mismatches){
                return 0;
            }
            for (int j=0; j<nmatches && number_of_variants+j < max_barcode_mismatches; j++){
                memcpy(corrected_codes+(number_of_variants+j)*barcode_code_length, code, barcode_code_length);
                corrected_codes[(number_of_variants+j)*barcode_code_length+i]=corrected_bases[j];
                indices[number_of_variants+j]=i;
            }
            number_of_variants+=nmatches;
        }
        if (number_of_variants > max_barcode_mismatches){
            return 0;
        }
    }
    return number_of_variants;
}
int find_variant_match(unsigned char *code, int sequence_index, unsigned char *corrected_bases ){
    const unsigned char index= sequence_index / 4;
    const unsigned char shift= 6-2*(sequence_index % 4);
    unsigned char mod_code[barcode_code_length];
    memcpy(mod_code, code, barcode_code_length);
    int number_of_variants=0; 
    for (unsigned char i=1; i<4; i++){
        unsigned char mod = code[index] ^ (i << shift );
        mod_code[index]=mod;
        if (g_hash_table_lookup(whitelist_hash, mod_code) != NULL){
            corrected_bases[number_of_variants++]=mod;
        }
    }
    return number_of_variants;
}
void process_pending_barcodes( data_structures *hashes, memory_pool_collection *pools, statistics *stats){
    unmatched_barcodes_features_block *current_entry_block=stats->unmatched_list.first_entry;
    while (current_entry_block != NULL){
        unsigned char *retcode=find_best_posterior_match(current_entry_block, number_of_features, MIN_POSTERIOR,stats);
        if (retcode != 0){
            unmatched_barcodes_features current_entry;
            read_unmatched_features_block(current_entry_block, &current_entry);
            char umi[umi_length+1];
            code2string(current_entry.umi,umi, umi_code_length);
            update_feature_counts_from_code(retcode, umi, current_entry.feature_index, hashes,pools);
        }
        current_entry_block=current_entry_block->next;
    }
}
unsigned char* find_best_posterior_match (unmatched_barcodes_features_block *entry_block, int number_of_features, double min_posterior,statistics *stats){
    unmatched_barcodes_features entry_struct;
    unmatched_barcodes_features *entry=&entry_struct;
    read_unmatched_features_block(entry_block, entry);
    if (!entry->feature_index){
        return 0;
    }
    double priors[max_barcode_mismatches+1];
    double evidence[max_barcode_mismatches+1];
    double total_evidence = 0.0;
    for (int i=0; i<entry->number_of_closest_barcodes; i++){
        //find counts for the barcode
        int total_counts=1;
        unsigned char* barcode =(entry->closest_barcodes) + i*barcode_code_length;
        feature_counts *counts=g_hash_table_lookup(whitelist_hash, barcode);
        if (counts != NULL){
            for (int j=1; j<=number_of_features; j++){
                total_counts+=counts->counts[j];
            }
            fprintf(stderr, "Total counts %d\n", total_counts);
        }
        priors[i]= pow(10,-0.1*(entry->Qscores[i]-33));
        evidence[i]=(double) total_counts * priors[i];
        total_evidence+=evidence[i];
    }
    for (int i = 0; i < entry->number_of_closest_barcodes; i++) {
        double posteriors = evidence[i] / total_evidence;
        if(posteriors > min_posterior){
            stats->pending_recovered++;
            return (entry->closest_barcodes) + i*barcode_code_length;
        }   
    }
    return 0;
}
int simpleCorrectFeature(char *line, feature_arrays *features, int maxN, int maxHammingDistance, int *hamming_distance){
    int best_feature=simple_search(features, line);
    *hamming_distance=0;
    if (!best_feature){
        const size_t length=strlen(line)-1;
        char buffer[(length+1) * (4 << ((maxN-1)*2))];
        char *corrected_seqs[ 4 << ((maxN-1)*2)];
        int nAlts=checkSequenceAndCorrectForN(line, corrected_seqs, buffer, length, maxN);
        if (nAlts == 1){
            return simple_hamming_search(features, corrected_seqs[0], maxHammingDistance,hamming_distance);
        }
        if (nAlts > 1){
            int bestHammingDistance=maxHammingDistance+1;
            for (int i=0; i<nAlts; i++){
                int feature_index=simple_search(features, corrected_seqs[i]);
                if (feature_index){
                    *hamming_distance=0;
                    return feature_index;
                }
            }
            for (int i=0; i<nAlts; i++){
                int myHammingDistance;
                int feature_index=simple_hamming_search(features, corrected_seqs[i], maxHammingDistance,&myHammingDistance);
                if (feature_index && myHammingDistance < bestHammingDistance){
                    bestHammingDistance=myHammingDistance;
                    best_feature=feature_index;
                }
                else if (feature_index && myHammingDistance == bestHammingDistance){
                    best_feature=0;
                }
            }

        }

    }
    return best_feature;
}
int checkAndCorrectFeature(char *line, feature_arrays *features,int maxHammingDistance, int nThreads, int *hamming_distance, char *matching_sequence, int maxN,char *ambiguous){
    const size_t length=strlen(line)-1;
    char buffer[(length+1) * (4 << ((maxN-1)*2))];
    char *corrected_seqs[ 4 << ((maxN-1)*2)];
    int hamming=0;
    //return ambiguous if the hamming distance is non-zero but the feature is zero
    //return ambiguous if there are too many Ns ie. nAlts is zero - distinguish this by setting hamming distance to maxHammingDistance+1
    int nAlts=checkSequenceAndCorrectForN(line, corrected_seqs, buffer, length, maxN);
    if (!nAlts){
        *hamming_distance=maxHammingDistance+1;
        *ambiguous=1;
        return 0;
    }
    if (nAlts == 1){
        char *myMatchingSequence;
        //no Ns and the barcode is good in terms of ACGT 
        int feature_index=find_feature_match_parallel(features, line, maxHammingDistance,nThreads,&hamming,&myMatchingSequence);
        *hamming_distance=hamming;
        if (feature_index){
            memcpy(matching_sequence, myMatchingSequence, features->feature_lengths[feature_index-1]);
            matching_sequence[features->feature_lengths[feature_index-1]]='\0';
        }
        else if (hamming <= maxHammingDistance){
            *ambiguous=1;
        }    
        return feature_index;
    }
    else{
        //if there are Ns in the barcode then we need to check all the possible sequences
        //and return if there is a unique match in the whitelist


        int number_of_matches=0;
        int myAmbiguous=0;
        int bestfeature_index=0;
        int bestHammingDistance=maxHammingDistance+1;
        char *bestMatchingSequence=0; 
        //int bestAlt=0;

        for (int i=0; i<nAlts; i++){
            char *myMatchingSequence;
            int feature_index=find_feature_match_parallel(features, corrected_seqs[i], maxHammingDistance,nThreads,&hamming,&myMatchingSequence);
            if (feature_index && !hamming){
                memcpy(matching_sequence, myMatchingSequence, features->feature_lengths[feature_index-1]);
                matching_sequence[features->feature_lengths[feature_index-1]]='\0';
                *hamming_distance=hamming;
                return feature_index;
            }
            else if (hamming < bestHammingDistance){
                if(number_of_matches > 1){
                    myAmbiguous=1;
                }
                else{
                    bestHammingDistance=hamming;
                    bestfeature_index=feature_index;
                    bestMatchingSequence=myMatchingSequence;
                    myAmbiguous=0;
                }
            }
            else if (hamming == bestHammingDistance){
                if( number_of_matches > 1 || feature_index != bestfeature_index){
                    myAmbiguous=1;
                }
            }
        }
        if (bestHammingDistance <= maxHammingDistance){
            *hamming_distance=bestHammingDistance;
            *ambiguous=myAmbiguous;
            if (!myAmbiguous && bestfeature_index){
                memcpy(matching_sequence, bestMatchingSequence, features->feature_lengths[bestfeature_index-1]);
                matching_sequence[features->feature_lengths[bestfeature_index-1]]='\0';
                return bestfeature_index;
            }
            return 0;
        }
        return 0;
    }
}
size_t barcode_code2number(unsigned char *code){
    unsigned char *key=(unsigned char*)g_hash_table_lookup(whitelist_hash, code);
    if (key == NULL){
        return 0;
    }
    return (key-whitelist)/barcode_code_length+1;
}
int checkAndCorrectBarcode(char **lines, int maxN, unsigned char feature_index, data_structures *hashes, memory_pool_collection *pools, statistics *stats){
    //fprintf(stderr, "Checking barcode %s\n", lines[1]);
    if (strlen(lines[1]) < barcode_length + umi_length){
        fprintf(stderr, "Error: Incomplete barcode %s\n", lines[1]);
        return 0;
    }
    if (!check_sequence(lines[1] + barcode_length, umi_length)){
        stats->nMismatches++;
        return 0;
    }
    //The return code should indicate whether we should calculate the barcode or not
    char buffer[(barcode_length + 1) * (4 << ((max_barcode_n-1)*2))];
    char *corrected_seqs[ 4 << ((max_barcode_n-1)*2)];
    char *candidateBarcode = lines[1];
    memset (buffer, 0, (barcode_length + 1) * (4 << ((max_barcode_n-1)*2)));
 
    //fprintf(stderr, "Checking barcode %s\n", candidateBarcode);    
    int nAlts=checkSequenceAndCorrectForN(candidateBarcode, corrected_seqs, buffer, barcode_length, maxN);
    if (!nAlts){
        stats->nMismatches++;
        return 0;
    }
    unsigned char code[barcode_code_length];
    if (nAlts == 1){
        //no Ns and the barcode is good in terms of ACGT 
        //check if the barcode is in the filtered whitelist 
        string2code(candidateBarcode, barcode_length, code); 
        if (g_hash_table_lookup(whitelist_hash, code)){
            update_feature_counts_from_code(code, lines[1]+barcode_length, feature_index, hashes, pools);
            stats->valid++;            
            return 1;
        }
        //if the barcode is not in the whitelist then try to correct it 
        unsigned char indices[max_barcode_mismatches+1];
        unsigned char corrected_codes[(max_barcode_mismatches+1)*barcode_code_length];
        stats->nMismatches++;
        //try to correct it by finding closest barcodes
        int nMatches=find_closest_barcodes(code,corrected_codes, indices);
        if (!nMatches){
            return 0;
        }

        //if unique closest barcode found then correct
        if (nMatches==1){
            char corrected_barcode[barcode_length+1];   
            code2string(corrected_codes, corrected_barcode, barcode_code_length);
            memcpy(lines[1], corrected_barcode, barcode_length);
            stats->recovered++;
            update_feature_counts_from_code(corrected_codes, lines[1]+barcode_length, feature_index, hashes, pools);
            stats->valid++;
            return 1;
        }
        else{
            //otherwise store it for later processing when we have counts to calculate the posterior probabilities
            if (nMatches > 1){
                unsigned char qscores[max_barcode_mismatches+1];
                for (int i=0; i<nMatches; i++){
                    qscores[i]=lines[3][indices[i]];
                }
                add_unmatched_barcode_store_feature(corrected_codes,lines[1]+barcode_length, qscores,feature_index,nMatches,pools,stats);
                stats->pending++;
            }
            return 0;
        }
        return 0;
    }
    else{
        stats->nMismatches++;
        //if there are Ns in the barcode then we need to check all the possible sequences
        //and return if there is a unique match in the whitelist
        int number_of_matches=0;
        unsigned char corrected_code[barcode_code_length];
        for (int i=0; i<nAlts; i++){
            unsigned char code[barcode_code_length];
            string2code(corrected_seqs[i], barcode_length, code);
            char corrected_barcode[barcode_length+1];   
            code2string(code, corrected_barcode, barcode_code_length);
            if (g_hash_table_lookup(whitelist_hash, code)){
                if (number_of_matches){
                    return 0;
                }
                else{
                    number_of_matches++;
                    memcpy(corrected_code, code, barcode_code_length);  
                }
            }
        }
        if (number_of_matches == 1){
            stats->recovered++;
            char corrected_barcode[barcode_length+1];   
            code2string(corrected_code, corrected_barcode, barcode_code_length);
            memcpy(lines[1], corrected_barcode, barcode_length);
            update_feature_counts_from_code(corrected_code, lines[1] + barcode_length,feature_index, hashes,pools);
            stats->valid++;
            return 1;
        }
        return 0;
    }
}
int find_number_of_R1_files(int positional_arg_count,char *fastqFilesR1string, char *fastqFilesR2string){
    int nFiles1=0;
    int nFiles2=0;
    if (positional_arg_count){
        return positional_arg_count/2;
    }
    for (int i=0; fastqFilesR1string[i]; i++){
        if (fastqFilesR1string[i] == ','){
            nFiles1++;
        }
    }
    for (int i=0; fastqFilesR2string[i]; i++){
        if (fastqFilesR2string[i] == ','){
            nFiles2++;
        }
    }
    if (nFiles1 != nFiles2){
        fprintf(stderr, "Error: Number of R1 and R2 files do not match\n");
        exit(EXIT_FAILURE);
    }
    return nFiles1+1;
}
int file_exists(const char *filename){
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}
const char* get_basename(const char* filename) {
    char *filepath = strdup(filename);
    //count the number of separators
    int number_of_separators=0;
    for (int i=0; filepath[i]; i++){
        if (filepath[i] == '/'){
            number_of_separators++;
        }
    }
    if (number_of_separators == 0){
        return filepath;
    }
    if (number_of_separators >= 1 && strlen(filepath) > 1){
        if (filepath[strlen(filepath)-1] == '/'){
            filepath[strlen(filepath)-1]='\0';
        }
    }
    const char* base = strrchr(filepath, '/');  // UNIX-like systems use '/'
    if (base) {
        return base + 1;  // Return everything after the last '/'
    }
    return filepath;  // If no '/' found, return the original string
}
// Function to return the basename up to the pattern "_L<numbers>_" or "_R1_"
char* extract_basename(const char* filepath, char *pattern) {
    //first get the basename
    const char* filename = get_basename(filepath);

    // try to find the "_R1_" pattern
    char *found_pattern = strstr(filename, pattern);
    if (found_pattern) {
        char *pointer=found_pattern-1;
        int digits_count=0;
        while (pointer > filename && isdigit(*pointer)){
            pointer--;
            digits_count++;
        }
        if (digits_count && *pointer=='L'){
            if (pointer > filename+1 && *(pointer-1) == '_'){
                found_pattern=pointer-1;
            }
            found_pattern=pointer-1;
        }
        size_t basename_len = found_pattern - filename;  // Calculate the length of the basename
        if (basename_len == 0) {
            return 0;
        }
        char *basename = (char*)malloc(basename_len + 1);  // Allocate memory for the basename
        if (!basename) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
        strncpy(basename, filename, basename_len);
        basename[basename_len] = '\0';  // Null-terminate the string
        return basename;
    }
    else{
        return 0;
    }
    
    // If neither pattern is found, return 0
    return 0;
}
void finalize_processing(feature_arrays *features, data_structures *hashes,  char *directory, int debug, memory_pool_collection *pools, statistics *stats){
    process_pending_barcodes(hashes, pools, stats);
    double elapsed_time = get_time_in_seconds() - stats->start_time;
    fprintf(stderr, "Finished processing %ld reads in %.2f seconds (%.1f thousand reads/second)\n", stats->number_of_reads, elapsed_time, stats->number_of_reads / (double)elapsed_time / 1000.0);

    //find size of whitelist_hash and non_whitelist_hash
    int total_feature_counts[features->number_of_features];
    int total_deduped_counts[features->number_of_features];
    int total_barcoded_counts[features->number_of_features];
    //coexpression matrix will store the total counts for each feature (also when there are no others - in the 0 field) in cells with feature i starting with 1
    //for simplicity we also hava a zero line so that all the indices are 1 based
    //store the total counts in the zero line
    int *coExpressionStorage = malloc((features->number_of_features + 1) * (features->number_of_features + 1) * sizeof(int));
    if (coExpressionStorage == NULL) {
        perror("Failed to allocate memory for coexpression matrix");
        exit(EXIT_FAILURE);
    }
    memset(coExpressionStorage, 0, (features->number_of_features + 1) * (features->number_of_features + 1) * sizeof(int));
    int **coExpression = malloc((features->number_of_features + 1) * sizeof(int *));
    if (coExpression == NULL) {
        perror("Failed to allocate memory for coexpression matrix");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < features->number_of_features + 1; i++) {
        coExpression[i] = coExpressionStorage + i * (features->number_of_features + 1);
    }
    int *histogramsStorage = malloc((features->number_of_features + 1) * (features->number_of_features + 1) * sizeof(int));
    if (histogramsStorage == NULL) {
        perror("Failed to allocate memory for coexpression histograms");
        exit(EXIT_FAILURE);
    }
    memset(histogramsStorage, 0, (features->number_of_features + 1) * (features->number_of_features + 1) * sizeof(int));
    int **coexpression_histograms = malloc((features->number_of_features + 1) * sizeof(int *));
    if (coexpression_histograms == NULL) {
        perror("Failed to allocate memory for coexpression histograms");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < features->number_of_features + 1; i++) {
        coexpression_histograms[i] = histogramsStorage + i * (features->number_of_features + 1);
    }
    //fprintf(stderr, "Number of reads matched to a feature %ld\n", valid);
    printFeatureCounts(features, total_deduped_counts, total_barcoded_counts, coExpression, coexpression_histograms, directory, hashes, stats);
    //fprintf(stderr, "Percentage of reads matched to a feature %.2f\n", (valid / (double)number_of_reads * 100.0));
    if (debug) {
        fprintf(stderr, "Number of mismatched reads that are matched to a nearest barcode unambiguously  %ld\n", stats->recovered);
        fprintf(stderr, "Number of reads pending priors before matching to a barcode %ld\n", stats->pending);
        fprintf(stderr, "Number of pending reads that were successfully matched to a barcode %ld\n", stats->pending_recovered);
    }

    print_feature_sequences(features, total_feature_counts, directory, hashes);
    char feature_stats_filename[FILENAME_LENGTH];
    sprintf(feature_stats_filename, "%s/features.txt", directory);
    FILE *feature_statsfp = fopen(feature_stats_filename, "w");
    if (feature_statsfp == NULL) {
        perror("Failed to open feature stats file");
        exit(EXIT_FAILURE);
    }
    fprintf(feature_statsfp, "Feature total_deduped_counts total_barcoded_counts total_feature_counts\n");

    char feature_coexpression_filename[FILENAME_LENGTH];
    sprintf(feature_coexpression_filename, "%s/feature_coexpression.txt", directory);
    FILE *feature_coexpression_fp = fopen(feature_coexpression_filename, "w");
    if (feature_coexpression_fp == NULL) {
        perror("Failed to open feature coexpression file");
        exit(EXIT_FAILURE);
    }
    int feature_printed[features->number_of_features];
    memset(feature_printed, 0, features->number_of_features * sizeof(int));
    int count = 0;
    for (int i = 0; i < features->number_of_features; i++) {
        if (feature_printed[i] == 0) {
            fprintf(feature_statsfp, "%s %d %d %d\n", features->feature_names[i], total_deduped_counts[count], total_barcoded_counts[count], total_feature_counts[count]);
            feature_printed[i] = 1;
            count++;
        }
    }
    for (int i = 0; i <= features->number_of_features; i++) {
        for (int j = 0; j < features->number_of_features + 1; j++) {
            fprintf(feature_coexpression_fp, "%d ", coExpression[i][j]);
        }
        fprintf(feature_coexpression_fp, "\n");
    }
    char feature_histograms_filename[FILENAME_LENGTH];
    sprintf(feature_histograms_filename, "%s/feature_histograms.txt", directory);
    FILE *feature_histograms_fp = fopen(feature_histograms_filename, "w");
    if (feature_histograms_fp == NULL) {
        perror("Failed to open feature histograms file");
        exit(EXIT_FAILURE);
    }
    for (int i = 1; i < features->number_of_features + 1; i++) {
        for (int j = 1; j < features->number_of_features + 1; j++) {
            fprintf(feature_histograms_fp, "%d ", coexpression_histograms[i][j]);
        }
        fprintf(feature_histograms_fp, "\n");
    }
    free(coExpression);
    free(coExpressionStorage);
    free(coexpression_histograms);
    free(histogramsStorage);
    fclose(feature_statsfp);
    fclose(feature_coexpression_fp);
    fclose(feature_histograms_fp);
}

void read_fastq_chunks(gzFile fastqFileR1gz, gzFile fastqFileR2gz, char *lines[8], char *buffer, char *done) {
    char *bufferptr = buffer;
    for (int i = 0; i < 8; i++) {
        gzFile gzfp = (i < 4) ? fastqFileR1gz : fastqFileR2gz;
        lines[i] = bufferptr;
        if (gzgets(gzfp, lines[i], LINE_LENGTH) != Z_NULL) {
            // Overwrite the newline character with null terminator
            const int length = strlen(lines[i]);
            if (lines[i][length - 1] != '\n') {
                fprintf(stderr, "Error: Incomplete record in the FASTQ file\n");
                exit(EXIT_FAILURE);
            }
            lines[i][length - 1] = '\0';
            bufferptr += length;
        } else {
            if (i == 0) {
                *done = 1;
                break;
            } else {
                fprintf(stderr, "Error: Incomplete record in the FASTQ file\n");
                exit(EXIT_FAILURE);
            }
        }
    }
}
void open_R1R2_files(const char *fastqFileR1, const char *fastqFileR2, gzFile *fastqFileR1gz, gzFile *fastqFileR2gz) {
    *fastqFileR1gz = gzopen(fastqFileR1, "r");
    if (*fastqFileR1gz == NULL) {
        fprintf(stderr, "Error: Unable to open R1 FASTQ file %s\n", fastqFileR1);
        exit(EXIT_FAILURE);
    }
    *fastqFileR2gz = gzopen(fastqFileR2, "r");
    if (*fastqFileR2gz == NULL) {
        fprintf(stderr, "Error: Unable to open R2 FASTQ file %s\n", fastqFileR2);
        gzclose(*fastqFileR1gz);
        exit(EXIT_FAILURE);
    }
}
void process_reads( data_structures *hashes, char *lines[8], feature_arrays *features, int maxHammingDistance, int nThreads, memory_pool_collection *pools, statistics *stats)
{
    stats->number_of_reads++;
    if (stats->number_of_reads % 1000000 == 0)
    {
        double elapsed_time = get_time_in_seconds() - stats->start_time;
        fprintf(stderr, "Processed %ld million reads in %.1f seconds\n", stats->number_of_reads / 1000000, elapsed_time);
    }
    // when directly acessing the features array - remember that the returned index is 1 based but the stored index is 0 based
    char matching_sequence[LINE_LENGTH];
    int hamming_distance = 0;
    unsigned char feature_index = 0;
    if (constant_offset >= 0)
    {
        feature_index = simpleCorrectFeature(lines[5] + constant_offset, features, max_feature_n, maxHammingDistance, &hamming_distance);
        if (feature_index)
        {
            // need to have feature_index -1 because the feature index is 1 based but array in struct is 0 based
            memcpy(matching_sequence, lines[5] + constant_offset, features->feature_lengths[feature_index - 1]);
            matching_sequence[features->feature_lengths[feature_index - 1]] = '\0';
        }
    }
    if (feature_index && hamming_distance)
    {
        char ambiguous = 0;
        const int myMaxHammingDistance = hamming_distance;
        int new_feature_index = checkAndCorrectFeature(lines[5], features, myMaxHammingDistance, nThreads, &hamming_distance, matching_sequence, max_feature_n, &ambiguous);
        if (hamming_distance < myMaxHammingDistance)
        {
            feature_index = new_feature_index; // even if zero keep the new feature index because that would mean ambiguity
        }
        // otherwise keep the old match even if there is an equal match outside of the constant region
    }
    else if (!feature_index)
    {
        char ambiguous = 0;
        feature_index = checkAndCorrectFeature(lines[5], features, maxHammingDistance, nThreads, &hamming_distance, matching_sequence, max_feature_n, &ambiguous);
    }
    if (feature_index)
    {
        // need to have feature_index -1 because the feature index is 1 based but array in struct is 0 based
        matching_sequence[features->feature_lengths[feature_index - 1]] = '\0';
        insert_feature_sequence(matching_sequence, feature_index, hamming_distance, hashes, pools);
        checkAndCorrectBarcode(lines, max_barcode_n, feature_index, hashes, pools, stats);
    }
    else
    {
        stats->nMismatches++;
        stats->total_unmatched_features++;
    }
}
void process_files_in_sample(char *directory, char **fastqFileR1, char **fastqFileR2, int nR1files, feature_arrays *features, int maxHammingDistance, int nThreads, int debug, memory_pool_collection *pools, statistics *stats, data_structures *hashes) {


    // Open the FASTQ files for reading
    gzFile fastqFileR1gz;
    gzFile fastqFileR2gz;

    mkdir_p(directory);
    for (int i = 0; i < nR1files; i++) {
        // Check if the output file directory/features_matrix.mtx exists
        open_R1R2_files(fastqFileR1[i], fastqFileR2[i], &fastqFileR1gz, &fastqFileR2gz);
        // Main loop for reading and processing the FASTQ files
        char done = 0;
        while (!done) {
            char buffer[LINE_LENGTH * 8];
            char *lines[8];
            read_fastq_chunks(fastqFileR1gz, fastqFileR2gz, lines, buffer, &done);
            // Read and process the second line from R1
            if (!done) {
                process_reads(hashes, lines, features, maxHammingDistance, nThreads, pools,stats);
            }
        }
        gzclose(fastqFileR1gz);
        gzclose(fastqFileR2gz);
        finalize_processing(features, hashes, directory, debug, pools,stats);
    }
}
void initialize_data_structures(data_structures *hashes){
    hashes->filtered_hash = g_hash_table_new(hash_int32, equal_int32);
    hashes->unique_features_match = g_hash_table_new(g_str_hash, g_str_equal);
    hashes->sequence_umi_hash = g_hash_table_new(hash_int64, equal_int64);
    hashes->neighbors_queue=malloc(sizeof(Queue));
    //check if the memory allocation was successful
    if (hashes->neighbors_queue == NULL) {
        perror("Failed to allocate memory for neighbors queue");
        exit(EXIT_FAILURE);
    }
    init_queue(hashes->neighbors_queue);
}
void destroy_data_structures(data_structures *hashes){
    g_hash_table_destroy(hashes->filtered_hash);
    g_hash_table_destroy(hashes->unique_features_match);
    g_hash_table_destroy(hashes->sequence_umi_hash);
    free_queue(hashes->neighbors_queue);
    free(hashes->neighbors_queue);
}

void find_sample_directory_name(char *directory, sample *sample, char *sample_directory, char sample_flag) {
    if (sample_flag) {
        char *name = extract_basename(sample->fastqFileR1[0], "_R1_");
        if (!name) {
            fprintf(stderr, "Error: Unable to extract basename from the FASTQ file name - expect _R1_ or _L<number>_R1_\n");
            exit(EXIT_FAILURE);
        }
        strcpy(sample_directory, directory);
        strcat(sample_directory, name);
        strcat(sample_directory, "/");
        free(name);
    } else {
        strcpy(sample_directory, directory);
    }
    fprintf(stderr, "Output directory %s\n", sample_directory);
}
void initialize_sample(sample *sample, int nFiles){
    fprintf(stderr, "Initializing sample with %d files\n", nFiles);
    char *block=malloc(2*nFiles*FILENAME_LENGTH);
    sample->fastqFileR1 = malloc(nFiles * sizeof(char *));
    for (int i=0; i<nFiles; i++){
        sample->fastqFileR1[i]=block+i*FILENAME_LENGTH;
        strcpy(sample->fastqFileR1[i], "hello");
    }
    sample->fastqFileR2 = malloc(nFiles * sizeof(char *));
    for (int i=0; i<nFiles; i++){
        sample->fastqFileR2[i]=block+(nFiles+i)*FILENAME_LENGTH;
        strcpy(sample->fastqFileR2[i], "hello");    
    } 
    sample->nFiles = 0;
}
void add_sample(sample *entry, char *fastqFileR1, char *fastqFileR2){
    int index = entry->nFiles;
    fprintf(stderr, "Adding sample %s %s to index %d\n", fastqFileR1, fastqFileR2, index);
    strcpy(entry->fastqFileR1[index], fastqFileR1);
    strcpy(entry->fastqFileR2[index], fastqFileR2);
    entry->nFiles++;
}
int organize_fastq_files_into_samples( char *fastqFileR1[], char *fastqFileR2[], int nR1files, sample **samples){
    int string_exists(char array[][FILENAME_LENGTH], int size, char *str) {
        for (int i = 0; i < size; i++) {
            if (strcmp(array[i], str) == 0) {
                return 1;  // String found
            }
        }
        return 0;  // String not found
    }
    int sample_size=0;
    char sample_names[nR1files][FILENAME_LENGTH];
    fprintf(stderr, "Number of R1 files %d\n", nR1files);
    for (int i=0; i<nR1files; i++){
        char *name1 = extract_basename(fastqFileR1[i],"_R1_");
        char *name2 = extract_basename(fastqFileR2[i],"_R2_");
        if (!name1 || !name2){
            fprintf(stderr, "Error: Unable to extract basename from the FASTQ file name  - expect _R1_ or _L<number>_R1_\n");
            exit(EXIT_FAILURE);
        }
        if (strcmp(name1, name2)){
            fprintf(stderr, "Error: R1 and R2 files do not match %s %s\n", name1, name2);
            exit(EXIT_FAILURE);
        }
        if (!string_exists(sample_names, sample_size, name1)){
            strcpy(sample_names[sample_size], name1);
            sample_size++;
        }
    }
    for (int i=0; i<sample_size; i++){
        for (int j=0; j<nR1files; j++){
            char *name1 = extract_basename(fastqFileR1[j],"_R1_");
            if (!strcmp(name1, sample_names[i])){
                add_sample(samples[i], fastqFileR1[j], fastqFileR2[j]);
            }
        }
    }
    return sample_size;
}
void get_fastq_filenames(int positional_arg_count, int argc, char *fastqFileR2[], char *argv[], char *fastqFileR1[], char *fastqFilesR1string, char *fastqFilesR2string){
    if (positional_arg_count > 0){
        for (int i = 0; optind < argc; i++, optind++){
            if (i % 2){
                fastqFileR2[i / 2] = argv[optind];
            }
            else{
                fastqFileR1[i / 2] = argv[optind];
            }
        }
    }
    else{
        // extract the FASTQ file names
        char *tokenR1 = strtok(fastqFilesR1string, ",");
        int fileCounts[2] = {0, 0};
        while (tokenR1)
        {
            fastqFileR1[fileCounts[0]] = tokenR1;
            fileCounts[0]++;
            tokenR1 = strtok(NULL, ",");
        }
        char *tokenR2 = strtok(fastqFilesR2string, ",");
        while (tokenR2)
        {
            fastqFileR2[fileCounts[1]] = tokenR2;
            fileCounts[1]++;
            tokenR2 = strtok(NULL, ",");
        }
        if (fileCounts[0] != fileCounts[1])
        {
            fprintf(stderr, "Error: Unequal number of R1 and R2 files\n");
            exit(EXIT_FAILURE);
        }
    }

}

int existing_output_skip(char keep_existing, char *directory){
    if (keep_existing && file_exists(directory)){
        char matrix_filename[4096];
        strcpy(matrix_filename, directory);
        strcat(matrix_filename, "features_matrix.mtx");
        if (file_exists(matrix_filename)){
            fprintf(stderr, "Matrix file %s found and skipping %s\n", matrix_filename, get_basename(directory));
            return 1;
        }
    }
    return 0;
}

void cleanup_sample(memory_pool_collection *pools, data_structures *hashes){        
    free_memory_pool_collection(pools);
    destroy_data_structures(hashes);
}
void calculate_sample_and_search_threads(int sample_size, int nThreads, int *sample_threads, int *search_threads){
    int search=(nThreads>=4)?4:1;
    int sample=nThreads/search;
    if (sample > sample_size){
        sample=sample_size;
    }
    *sample_threads=sample;
    *search_threads=search;
}

int main(int argc, char *argv[])
{
    omp_set_nested(1);
    initseq2Code();
    initcode2seq();
    initdiff2hamming(diff2Hamming);
    feature_arrays *features=0;
    char *fastqFilesR1string=0;
    char *fastqFilesR2string=0;
    int nThreads=4;
    int maxHammingDistance=1;
    char directory[LINE_LENGTH];
    strcpy (directory, "features/");
    int c=0;
    char whitelist_filename[4096];
    char sample_flag=1;
    char keep_existing=0;   
    static struct option long_options[] = {
        {"whitelist", required_argument, 0, 'w'},
        {"featurelist", required_argument, 0, 'f'},
        {"maxHammingDistance", required_argument, 0, 'm'},
        {"constant_offset", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"stringency", required_argument, 0, 's'},
        {"min_counts", required_argument, 0, 'i'},
        {"umi_length",required_argument , 0, 'u'},
        {"barcode_length",required_argument , 0, 'b'},
        {"directory", required_argument, 0, 'd'},
        {"debug", no_argument, 0, 'v'},
        {"as_named", no_argument, 0, 'a'},
        {"keep_existing", no_argument, 0, 'k'},
        {"featureR1", required_argument, 0, 0},
        {"featureR2", required_argument, 0, 1},
        {"max_barcode_mismatches", required_argument, 0, 2},
        {"feature_n", required_argument, 0, 3},
        {"barcode_n", required_argument, 0, 4},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    while ((c = getopt_long(argc, argv, "w:b:f:m:s:i:t:o:d:u:vak", long_options, &option_index)) != -1) {
        switch (c) {
            case 'w':
                strcpy(whitelist_filename, optarg);
                break;
            case 'b':
                barcode_length=atoi(optarg);
                barcode_code_length=(barcode_length+3)/4;
                break;
            case 'f':
                features=readTagsFile(optarg);
                number_of_features=features->number_of_features;
                maximum_feature_length=features->max_length;
                feature_code_length=(maximum_feature_length+3)/4;
                break;
            case 'm':
                maxHammingDistance=atoi(optarg);
                break;
            case 's':
                stringency=atoi(optarg);
                break;
            case 'i':
                min_counts=atoi(optarg);
                break;
            case 'd':
                strcpy(directory, optarg);
                //check if directory ends in /
                if (directory[strlen(directory)-1] != '/'){
                    strcat(directory, "/");
                }
                break;
            case 'o':
                constant_offset=atoi(optarg);
                break;
            case 't':
                nThreads=atoi(optarg);
                break;
            case 'u':
                umi_length=atoi(optarg);
                umi_code_length=(umi_length+3)/4;
                break;
            case 'v':  // Handle the -v flag for enabling debug mode
                debug = 1;
                break;
            case 'a':
                sample_flag=0;
                break;
            case 'k':
                keep_existing=1;
                break;
            case 0:
                fastqFilesR1string = malloc(strlen(optarg) + 1);
                if (fastqFilesR1string == NULL) {
                    fprintf(stderr, "Error: Unable to allocate memory for fastqFilesR1string\n");
                    exit(EXIT_FAILURE);
                }
                strcpy(fastqFilesR1string, optarg);
                fprintf(stderr, "Read fastqFilesR1String %s %s\n", optarg, fastqFilesR1string);
                break;
            case 1:
                fastqFilesR2string = malloc(strlen(optarg) + 1);
                if (fastqFilesR2string == NULL) {
                    fprintf(stderr, "Error: Unable to allocate memory for fastqFilesR2string\n");
                    exit(EXIT_FAILURE);
                }
                strcpy(fastqFilesR2string, optarg);
                break;
            case 2:
                max_barcode_mismatches=atoi(optarg);
                break;    
            case 3:
                max_feature_n=atoi(optarg);
                break;
            case 4:
                max_barcode_n=atoi(optarg);
                break;
            default:
                // print usage
                fprintf(stderr, "Usage: %s -w cell_barcodes_whitelist_file -f feature_barcode_files -s stringency -i stringency_min_counts -d output_directory -o constant_offset -t number_of_threads -u umi_length -v (verbose output) -a (do not put in directory by sample name) --featureR1 <comma separate list of R1 fastq files> --featureR2  <comma separate list of R1 fastq files> --barcode_N <max Ns in cell barcodes> --feature_N <max Ns in feature barcodes> --max_barcode_mismatches <max mismatch barcode candidates>  <feature files as arguments if --featureR1 not used>\n", argv[0]);
                return 1;
        }
    }

    int positional_arg_count = argc - optind;
    //sanity check the parameters
    if (positional_arg_count == 0 && (!fastqFilesR1string || !fastqFilesR2string ) ){
        fprintf(stderr, "Error: No FASTQ files provided\n");
        exit(EXIT_FAILURE);
    }
    if (positional_arg_count > 0) {
        if (fastqFilesR1string || fastqFilesR2string) {
            fprintf(stderr, "Error: FASTQ files provided both as positional arguments and as options\n");
            exit(EXIT_FAILURE);
        }
        if (positional_arg_count % 2 != 0) {
            fprintf(stderr, "Error: Unequal number of R1 and R2 files\n");
            exit(EXIT_FAILURE);
        }
    }
    int nR1files=find_number_of_R1_files(positional_arg_count,fastqFilesR1string, fastqFilesR2string);
    char *fastqFileR1[nR1files];
    char *fastqFileR2[nR1files];
    whitelist_hash = g_hash_table_new(hash_int32, equal_int32 );
    readWhiteList(whitelist_filename, whitelist_hash);
    initialize_unit_sizes();
    get_fastq_filenames(positional_arg_count, argc, fastqFileR2, argv, fastqFileR1, fastqFilesR1string, fastqFilesR2string);
    sample *samples[nR1files];
    sample *block=0;
    int nSamples=0; 
    if (sample_flag){
        block = malloc(sizeof(sample)*nR1files);
        memset(block, 0, sizeof(sample)*nR1files);
        for (int i=0; i<nR1files; i++){
            samples[i]=block+i;
            initialize_sample(samples[i], nR1files);
        }
        nSamples=organize_fastq_files_into_samples(fastqFileR1, fastqFileR2, nR1files, samples);
    }
    else{
        nSamples=1;
        initialize_sample(samples[0], nR1files);
        for (int i=0; i<nR1files; i++){
            add_sample(samples[0], fastqFileR1[i], fastqFileR2[i]);
        }
    }
    fprintf(stderr, "Number of samples %d\n", nSamples); 
    int sample_threads=1;
    int search_threads=4;
    calculate_sample_and_search_threads(nSamples, nThreads, &sample_threads, &search_threads);
    fprintf(stderr, "Sample threads %d Search threads %d\n", sample_threads, search_threads);
    if (nSamples > 1 && sample_threads > 1){   
        #pragma omp parallel for num_threads(sample_threads)  
        for (int i=0; i<nSamples; i++){
            char sample_directory[FILENAME_LENGTH];
            memory_pool_collection *pools=initialize_memory_pool_collection();
            statistics stats;
            data_structures hashes;
            initialize_data_structures(&hashes);
            initialize_statistics(&stats);
            find_sample_directory_name(directory, samples[i], sample_directory, sample_flag);
            fprintf(stderr, "Processing sample directory %s\n", sample_directory);
            if (existing_output_skip(keep_existing, sample_directory)) continue;
            process_files_in_sample(sample_directory, samples[i]->fastqFileR1, samples[i]->fastqFileR2, samples[i]->nFiles, features, maxHammingDistance, search_threads, debug, pools, &stats,&hashes);
            cleanup_sample(pools, &hashes);
        }
    }
    else{
         for (int i=0; i<nSamples; i++){
            char sample_directory[FILENAME_LENGTH];
            memory_pool_collection *pools=initialize_memory_pool_collection();
            statistics stats;
            data_structures hashes;
            initialize_data_structures(&hashes);
            initialize_statistics(&stats);
            find_sample_directory_name(directory, samples[i], sample_directory, sample_flag);
            fprintf(stderr, "Processing sample %s\n", get_basename(sample_directory));
            if (existing_output_skip(keep_existing, sample_directory)) continue;
            process_files_in_sample(sample_directory, samples[i]->fastqFileR1, samples[i]->fastqFileR2, samples[i]->nFiles, features, maxHammingDistance, search_threads, debug, pools, &stats,&hashes);
            cleanup_sample(pools, &hashes);
        }       
    }
    g_hash_table_destroy(whitelist_hash);
    free_feature_arrays(features);
    return 0;
}

