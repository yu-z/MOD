#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

# take the list of commands as input and run through qsub
############# qsub parameter ################
$PM =
'#!/bin/sh
#$ -cwd 
#$ -o /u/scratch2/j/jxzhai/qsub/job-output-files/
#$ -M jxzhai@ucla.edu
#$ -m a
#$ -V
#$ -S /bin/bash
#$ -l h_data=2G,h_rt=23:55:00
#$ -e /u/scratch2/j/jxzhai/qsub/job-error-files/
#$ -pe shared 4
';

############# BSMAP parameter ########################
$refG = '/u/home/j/jxzhai/MCDB_folder/ALLDATA/Genomic/TAIR10/GenomicSEQ/Ath_ChrAll.fa';
#single_end: using 8 threads [-p 8], and allow up to 1 multiple hits [-w 1], map to two forward possible strands [-n 0], allow 2 mismatches [-v 2], 
# -u report unmapped reads
#-g  <int>   gap size, BSMAP only allow 1 continuous gap (insertion or deletion) with up to 3 nucleotides. default=0
            # gaps will not be allowed within 6nt of the read edges.
            # the number of mismatches of gapped algnment is calculated as #gap_size+#mismatches+1
$BSMAP_PM = '-p 4 -w 1 -n 0 -v 0.04 -u -g 0 '; # allow 4% mismatch

#################### open folder contains all fastq files#################
my $USAGE = "\nUSAGE: multiple_BSMAP.pl 
                                   -seqdir directory_with_fastq_files
                                   ";
my $options = {};
GetOptions($options, "-seqdir=s"); #, "-out=s" 
die $USAGE unless defined ($options->{seqdir});

############################# Grobal Variables #############################
my $seqdir = $options->{seqdir};

opendir (DIR, $seqdir) or die "Couldn't open $seqdir: $!\n";

my $command;
while (defined(my $seqfile = readdir(DIR))) {
    unless ($seqfile =~ /(.*)\.fastq\.gz$/) {next;}
    $lib_name = $1;
    if ($seqfile =~ /\.gz$/) {
    open(SEQ1, "gunzip -c $seqdir$seqfile |") || die "can't open pipe to $seqfile";
    }
    else {
    open(SEQ1, "< $seqdir$seqfile") || die "can't open $seqfile";
    }

$command =
"
bsmap -a $seqdir$seqfile -d $refG -o $lib_name\.bam $BSMAP_PM

python /u/home/j/jxzhai/install/bsmap-2.74/methratio_alt.py --Steve --ref=$refG --out=$lib_name -u -z -r $lib_name\.bam

# perl ~/Script/Analysis/BS-seq/bsmap2wiggle_v1.0.pl $lib_name.gz

# gzip $lib_name*.wig

perl ~/Script/Analysis/BS-seq/BSmap2cytosine.pl --input_file $lib_name.gz --reference_cytosine ~/MCDB_folder/ALLDATA/Genomic/TAIR10/C_position/TAIR10_v2.cytosine.gz

perl ~/Script/Analysis/BS-seq/Cytosine_2_100bp_v2.pl $lib_name.cytosine.gz

rm $lib_name.cytosine.gz
rm $lib_name\.bam
rm $lib_name\.bam\.bai
rm $lib_name\.gz

";
    print "$command\n";
    open(SH, ">SE\_$lib_name\_BSMAP.sh");
    print SH "$PM";
    print SH "$command\n";
    system "qsub SE\_$lib_name\_BSMAP.sh\n";
    close SH;
    system "rm SE\_$lib_name\_BSMAP.sh";
}