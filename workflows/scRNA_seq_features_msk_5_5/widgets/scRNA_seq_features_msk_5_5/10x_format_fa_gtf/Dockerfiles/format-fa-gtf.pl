#!/usr/bin/perl
use File::Basename;
my $fa="$ENV{INPUTFA}";
my $outputfa="$ENV{OUTPUTFA}";
my $gtf="$ENV{INPUTGTF}";
my $outputgtf="$ENV{OUTPUTGTF}";
my $overwrite="$ENV{overwrite}";
my $featuregtf="$ENV{INPUTFGTF}";
my $featurefa="$ENV{INPUTFFA}";


if (-e $outputfa && !$overwrite){
    print STDERR "$outputfa already exists and overwrite is not enable - skipping fa generation\n";
}
else{
    print STDERR "working on $fa\n";
    makeDirectory($outputfa);
    modifyChrNames($fa,$outputfa);
    if ($featurefa){
        system("cat $featurefa >> $outputfa");
    }
}
if (-e $outputgtf && !$overwrite){
    print STDERR "$outputgtf already exists and overwrite is not enable - skipping gtf generation\n";
}
else{
    print STDERR "working on $gtf\n";
    makeDirectory($outputgtf);
    %passedFilter=generateFilterForGTF($gtf);
    modifyGTF($gtf,$outputgtf);
    if ($featuregtf){
        system("cat $featuregtf >> $outputgtf");
    }
}

sub makeDirectory{
    my($path)=@_;
    my $directory = dirname( $path );
    unless(-d "$directory"){
        system("mkdir -p $directory")
    }

}

sub modifyGTF{
    my($gtf,$outputgtf)=@_;
    open (my $fp,"$gtf") || die "can't open $gtf";
    open (my $outfp,">$outputgtf") || die "can't open $outputgtf";
    while (defined(my $line=<$fp>)){
        if(substr($line,0,1) eq "#"){
            print $outfp  "$line"
        }
        else{
            chomp($line);
            my (@parts)=split(/\t/,$line);
            #check if Y chromosome and PAR
            if ($parts[0] eq "chrY" && $parts[8] =~ /tag \"PAR\"/){next;}
            my (@infoParts)=split(/\;/,$parts[8]);
            if($infoParts[0] =~ /(ENSG[0-9]+)\.([0-9]+)/){
                my $geneID = $1;
                if ($passedFilter{$geneID}){
                    $infoParts[0] = "gene_id \"$geneID\"; gene_version \"$2\"";
                    if($infoParts[1] =~ /(ENST[0-9]+)\.([0-9]+)/){
                      $infoParts[1] = " transcript_id \"$1\"; transcript_version \"$2\"";
                    }
                    $parts[8]=join(";",@infoParts);
                    $parts[8] =~ s/(exon_id \"ENSE[0-9]+)\.([0-9]+\")/$1\"\; exon_version \"$2/;
                    printf $outfp "%s;\n",join("\t",@parts);

                }
 
            }        
        }
    }    

}
sub generateFilterForGTF{
    my($gtf)=@_;
    my %goodGenes;
    my $biotypes="(protein_coding|lncRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene)";    
    open (my $fp,"$gtf") || die "can't open $gtf";
    while (defined(my $line=<$fp>)){
        chomp($line);
        my (@parts)=split(/\t/,$line);
        if($parts[2] ne "transcript"){next}
        my (@infoParts)=split(/\;/,$parts[8]);
        if($infoParts[0] =~ /(ENSG[0-9]+)\./){
            my $geneID = $1;
            if ($goodGenes{$geneID}){next}
            if ( $parts[8] =~ /gene_type \"$biotypes\"/ && $parts[8] =~ /transcript_type \"$biotypes\"/ && $parts[8] !~ /tag \"PAR\"/ && $parts[8] !~ /tag \"readthrough_transcript\"/ ) {
                $goodGenes{$geneID}++;
            } 
        }
    }
    return %goodGenes;
}

sub modifyChrNames{
    my($fa,$outputfa)=@_;
    open (my $fp,"$fa") || die "can't open $fa";
    open (my $outfp,">$outputfa") || die "can't open $outputfa";
    while (defined(my $line=<$fp>)){
        if(substr($line,0,1) eq '>'){  
            my ($chr)=split(' ',$line);
            if($chr =~ /^\>([0-9XY]+)/){
                print $outfp ">chr$1 $1\n";
            }
            elsif(substr($chr,1,2) eq "MT"){
                print $outfp ">chrM MT\n";
            }
            else{
              printf $outfp "$chr %s\n",substr($chr,1);
            }
        }
        else{
            print $outfp "$line";
        }
    }

}