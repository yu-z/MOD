#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

# take the list of SRR as input and parse the BSseeker2 log file ".bs_seeker2_log"

#################### open folder contains all fastq files#################
my $USAGE = "\nUSAGE: Parse_BSMAP_log.pl 
                                   -seqdir directory_with_BSseeker2_output_files
                                   -srr file_that_contain_list_of_SRR
                                   ";
my $options = {};
GetOptions($options, "-seqdir=s", "-srr=s" ); #, "-out=s" 
die $USAGE unless defined ($options->{seqdir});
# die $USAGE unless defined ($options->{srr});

############################# Grobal Variables #############################
my $seqdir = $options->{seqdir};
my $srr = $options->{srr};

opendir (DIR, $seqdir) or die "Couldn't open $seqdir: $!\n";
%hash;

while (defined(my $seqfile = readdir(DIR))) {
    unless ($seqfile =~ /(.*)\.CG.100.gz/) {next;}
    $lib_name = $1;
    my @line = `head $seqdir$lib_name.BSMAP_log`;
    unless ($line[0] =~ m/total (\d+) valid mappings, (\d+) reads tossed by Steve filter, (\d+) covered cytosines, average coverage: (.+) fold./) {next;}
    my $uniquely_mapped_reads = $1;
    my $steve = $2;
    my $coverage = $4;
    my %m_ratio_Chr1_5;
    my %m_ratio_Chloroplast;
    my %passed_bin;
    
    foreach my $type ("CG","CHG","CHH") {
        my $file = "$seqdir$lib_name.$type.100";
        open($type, "gunzip -c $file |") || die "can't open pipe to $file";
        my %mC_count;
        my %uC_count;
        my $total_bin = 0;
        while (<$type>){
            chomp;
            next if (m/chr/i);
            @data = split/\t/;
            $mC_count{$data[0]} += $data[2];
            $uC_count{$data[0]} += $data[3];
            $passed_bin{$type} ++ if ($data[5]>=4);
            $total_bin ++;
        }
        my $mC_genome = $mC_count{"1"} + $mC_count{"2"} + $mC_count{"3"} + $mC_count{"4"} + $mC_count{"5"};
        my $uC_genome = $uC_count{"1"} + $uC_count{"2"} + $uC_count{"3"} + $uC_count{"4"} + $uC_count{"5"};
        my $mC_chloroplast = $mC_count{"C"};
        my $uC_chloroplast = $uC_count{"C"};
        $m_ratio_Chr1_5{$type} = sprintf("%.5f", $mC_genome/($uC_genome + $mC_genome));
        # print "$uC_genome\t$mC_genome\n";
        $m_ratio_Chloroplast{$type} = sprintf("%.5f", $mC_chloroplast/($mC_chloroplast + $uC_chloroplast)); 
        # print "\t$type\t$m_ratio\t$m_ratio_C\t$passed_bin\t$total_bin";
    }
    $hash{$lib_name} = "$lib_name\t$uniquely_mapped_reads\t$steve\t$coverage\t$m_ratio_Chloroplast{\"CG\"}\t$m_ratio_Chloroplast{\"CHG\"}\t$m_ratio_Chloroplast{\"CHH\"}\t$m_ratio_Chr1_5{\"CG\"}\t$m_ratio_Chr1_5{\"CHG\"}\t$m_ratio_Chr1_5{\"CHH\"}\t$passed_bin{\"CG\"}\t$passed_bin{\"CHG\"}\t$passed_bin{\"CHH\"}\t$total_bin";
    print "$lib_name\t$uniquely_mapped_reads\t$steve\t$coverage\t$m_ratio_Chloroplast{\"CG\"}\t$m_ratio_Chloroplast{\"CHG\"}\t$m_ratio_Chloroplast{\"CHH\"}\t$m_ratio_Chr1_5{\"CG\"}\t$m_ratio_Chr1_5{\"CHG\"}\t$m_ratio_Chr1_5{\"CHH\"}\t$passed_bin{\"CG\"}\t$passed_bin{\"CHG\"}\t$passed_bin{\"CHH\"}\t$total_bin\n"; 
    undef %m_ratio_Chr1_5;
    undef %m_ratio_Chloroplast;
    undef %passed_bin;
}

open(IN, $srr) or die "Couldn't open $srr: $!\n";;
open OUT, "summary.txt";

while (<IN>) {
    chomp;
    print OUT "$_\t$hash{$_}\n";

}

