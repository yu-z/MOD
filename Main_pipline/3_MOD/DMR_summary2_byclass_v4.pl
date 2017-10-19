#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

#################### open folder contains all input files#################
my $USAGE = "\nUSAGE: DMR_summary2.pl 
                                   -seqdir DMR_folder
                                   -wtlist file containing the list of WT library IDs
                                   -mulist file containing the list of muant library IDs
                                   -type CG,CHG, or CHH
                                   ";
my $options = {};
GetOptions($options, "-seqdir=s", "-wtlist=s", "-mulist=s", "-type=s"); #, "-out=s" 
die $USAGE unless defined ($options->{seqdir});
die $USAGE unless defined ($options->{wtlist});
die $USAGE unless defined ($options->{mulist});
die $USAGE unless defined ($options->{type});

############################# Grobal Variables #############################
my $seqdir = $options->{seqdir};
my $wtlist = $options->{wtlist};
my $mulist = $options->{mulist};
my $run_type = $options->{type};

opendir (DIR, $seqdir) or die "Couldn't open $seqdir: $!\n";


# obtain the list of WT libraries and save as hash
open (WT,$wtlist);
my %WT_hash;
my @WT_list;
# my %mutant_hash;
while (<WT>){
    chomp;
    next if (length($_)<1);
    $WT_hash{$_}++;
    push (@WT_list, $_);
}

open (MU,$mulist);
my %MU_hash;
my @MU_list;

while (<MU>){
    chomp;
    next if (length($_)<1);
    $MU_hash{$_}++;
    push (@MU_list, $_);
}

$WT_size = keys %WT_hash;
$MU_size = keys %MU_hash;

# @cutoff_list = (1, 5, 10, 15, 20, 25);
@cutoff_list = (1..54);

$lib_cutoff = $MU_size;

print "total number of WT libraries: $WT_size , total number of MU libraries: $MU_size , \n";


my %mutant_hash;
my %overlap_hash;

while (defined(my $seqfile = readdir(DIR))) {
    next if (($seqfile =~ m/^\.$/) or ($seqfile =~ m/^\.\.$/));
    next if !($seqfile =~ m/DMR/);
    next if !($seqfile =~ m/^$run_type\./);

    
    $seqfile =~ /(.*)\.DMR\.(.*)_VS_(.*)\.txt/;
    my $type = $1;
    my $lib1 = $2; $lib1 =~ s/#\d+//;
    my $lib2 = $3; $lib2 =~ s/#\d+//;
    
    if (!(exists $MU_hash{$lib1}) && !(exists $MU_hash{$lib2})) {next;}
    if ((exists $WT_hash{$lib1}) && (exists $WT_hash{$lib2})) {
        # if (exists $WT_hash{$lib2}) {next;}
        open(FH,"$seqdir$seqfile") or die "cannot open file $!";
        while (<FH>) {
            chomp;
            next if (m/chr/i);
            my $class;
            my @data = split /\t/;
            my $bin_key = join ('_',$data[0],$data[1]);
            my($c1,$ct1,$c2,$ct2) = ($data[2],$data[3],$data[4],$data[5]);
            if (($c1/($c1+$ct1)) > ($c2/($c2+$ct2))) {
                $mutant_hash{$lib2}{$lib1}{$type}{"hypo"}{$bin_key}++; 
                $mutant_hash{$lib1}{$lib2}{$type}{"hyper"}{$bin_key}++; 

                $overlap_hash{$lib2}{$type}{"hypo"}{$bin_key}++;  
                $overlap_hash{$lib1}{$type}{"hyper"}{$bin_key}++;  
            }
            else {
                $mutant_hash{$lib2}{$lib1}{$type}{"hyper"}{$bin_key}++; 
                $mutant_hash{$lib1}{$lib2}{$type}{"hypo"}{$bin_key}++;

                $overlap_hash{$lib2}{$type}{"hyper"}{$bin_key}++;  
                $overlap_hash{$lib1}{$type}{"hypo"}{$bin_key}++;  
            }
            ###### count of bin_key occurrence
            # $overlap_hash{$mutant}{$type}{$class}{$bin_key}++;          
        }
        close FH;
        next;
    }
    elsif ((exists $WT_hash{$lib1}) && !(exists $WT_hash{$lib2})) {
        # if (exists $WT_hash{$lib2}) {next;}
        next if !(exists $MU_hash{$lib2});
        $WT = $lib1; $mutant = $lib2;
    }
    elsif (!(exists $WT_hash{$lib1}) && (exists $WT_hash{$lib2})) {
        # if (exists $WT_hash{$lib1}) {next;}
        next if !(exists $MU_hash{$lib1});
        $WT = $lib2; $mutant = $lib1;
    }
    else {
        next;
    }
    # print "$WT\t$mutant\n";

    open(FH,"$seqdir$seqfile") or die "cannot open file $!";
    while (<FH>) {
        chomp;
        next if (m/chr/i);
        my $class;
        my @data = split /\t/;
        my $bin_key = join ('_',$data[0],$data[1]);
        my($c1,$ct1,$c2,$ct2) = ($data[2],$data[3],$data[4],$data[5]);
        if (($c1/($c1+$ct1)) > ($c2/($c2+$ct2))) {
            if ($WT eq $lib1 ) {$class = "hypo" ;} 
            if ($WT eq $lib2 ) {$class = "hyper" ;} 
        }
        else {
            if ($WT eq $lib1) {$class = "hyper" ;} 
            if ($WT eq $lib2) {$class = "hypo" ;} 
        }
        ###### count of bin_key occurrence
        $mutant_hash{$mutant}{$WT}{$type}{$class}{$bin_key}++;
        $overlap_hash{$mutant}{$type}{$class}{$bin_key}++;        
    }
    close FH;
 
}


########### separate hyper and hypo DMRs #########


# my %overlap;
# my %total_N_hash;
# my %total_N_hash_filter;
# my %total_N;
# my %total_Nfilter;



foreach my $class ("hypo","hyper")
{
    # foreach my $type ("CG","CHG","CHH")
    my $type = $run_type;
    {
        open($type, ">$class\_$type\_DMR_summary.txt");

        print $type "class\tcontext\tmutant_lib";
        foreach my $cutoff (@cutoff_list){
            print $type "\tintersect\_$cutoff";    
        }
        foreach my $lib_wt (@WT_list){
            print $type "\t$lib_wt";
        }
        print $type "\n";

        foreach my $lib_mu (@MU_list){
            print $type "$class\t$type\t$lib_mu";

            # count of DMR that passes cutoff
            my %intersect;
            foreach my $bin (keys %{$overlap_hash{$lib_mu}{$type}{$class}}){
                foreach my $cutoff (@cutoff_list) {
                    if ($overlap_hash{$lib_mu}{$type}{$class}{$bin} >= $cutoff) {
                    $intersect{$cutoff}++;    
                    }
                }
            }
            foreach my $cutoff (@cutoff_list){
                my $value = $intersect{$cutoff} + 0;
                print $type "\t$value";    
            }
            undef %intersect;

            # count of DMR for each WT
            foreach my $lib_wt (@WT_list){
                my $count_lib = keys %{$mutant_hash{$lib_mu}{$lib_wt}{$type}{$class}};
                print  $type "\t$count_lib";
            }
            print $type "\n";
        }
    # die;
    }
}
