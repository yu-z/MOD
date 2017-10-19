#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;
use Statistics::R;

#################### open folder contains all input files#################
my $USAGE = "\nUSAGE: DMR_summary2.pl 
                                   -seqdir DMR_folder
                                   -wtlist file containing the list of WT library IDs
                                   -mulist file containing the list of muant library IDs
                                   ";
my $options = {};
GetOptions($options, "-seqdir=s", "-wtlist=s", "-mulist=s"); #, "-out=s" 
die $USAGE unless defined ($options->{seqdir});
die $USAGE unless defined ($options->{wtlist});
die $USAGE unless defined ($options->{mulist});

############################# Grobal Variables #############################
my $seqdir = $options->{seqdir};
my $wtlist = $options->{wtlist};
my $mulist = $options->{mulist};
opendir (DIR, $seqdir) or die "Couldn't open $seqdir: $!\n";

# obtain the list of WT libraries and save as hash
open (WT,$wtlist);
my %WT_hash;
while (<WT>){
    chomp;
    next if (length($_)<1);
    $WT_hash{$_}++;
}

open (MU,$mulist);
my %MU_hash;

while (<MU>){
    chomp;
    next if (length($_)<1);
    $MU_hash{$_}++;
}

$WT_size = keys %WT_hash;
$MU_size = keys %MU_hash;

$cutoff = 33;

$min_sample_size = 0;

$lib_cutoff = $MU_size;

print "total number of WT libraries: $WT_size , cutoff set at $cutoff , total number of MU libraries: $MU_size , \n";

my %mutant_hash;

while (defined(my $seqfile = readdir(DIR))) {
    next if (($seqfile =~ m/^\.$/) or ($seqfile =~ m/^\.\.$/));
    next if !($seqfile =~ m/DMR/);
    
    $seqfile =~ /(.*)\.DMR\.(.*)_VS_(.*)\.txt/;
    my $type = $1;
    my $lib1 = $2; $lib1 =~ s/#\d+//;
    my $lib2 = $3; $lib2 =~ s/#\d+//;
    
    if (exists $WT_hash{$lib1}) {
        if (exists $WT_hash{$lib2}) {next;}
        next if !(exists $MU_hash{$lib2});
        $WT = $lib1; $mutant = $lib2;
    }
    elsif (exists $WT_hash{$lib2}) {
        if (exists $WT_hash{$lib1}) {next;}
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
            if (exists $WT_hash{$lib1}) {$class = "hypo" ;} 
            if (exists $WT_hash{$lib2}) {$class = "hyper" ;} 
        }
        else {
            if (exists $WT_hash{$lib1}) {$class = "hyper" ;} 
            if (exists $WT_hash{$lib2}) {$class = "hypo" ;} 
        }
        ###### count of bin_key occurrence
        $mutant_hash{$mutant}{$type}{$class}{$bin_key}++;        
    }
    close FH;
 
}

########### separate hyper and hypo DMRs #########

my @array;
my $n=0;


open (REF,"/u/home/j/jxzhai/MCDB_folder/Projects/BigdataBSseq/Round2_DMR/ref.lst");

while (<REF>){
    chomp;
    ($ID,$name) = split /\t/;
    $hash_ID{$ID} = $name;
    # print "$ID\t$name\n";
}


foreach my $lib1 (sort keys %mutant_hash){
    $n++;
    $array[$n][0] = "$lib1";
    $array[0][$n] = "$lib1";
}
print "total number of mutant libraries: $n\n";


my %intersect;
my %overlap;
my %total_N_hash;
my %total_N_hash_filter;
my %total_N;
my %total_Nfilter;


foreach my $class ("hypo","hyper")
{
    foreach my $type ("CG","CHG","CHH")
    {
        foreach my $lib1 (sort keys %mutant_hash){
            foreach my $bin (keys %{$mutant_hash{$lib1}{$type}{$class}}){
                if ($mutant_hash{$lib1}{$type}{$class}{$bin} >= $cutoff) {
                $intersect{$lib1}{$type}{$class}{$bin} = $mutant_hash{$lib1}{$type}{$class}{$bin};
                $total_N_hash_filter{$type}{$class}{$bin}++;
                }
                $total_N_hash{$type}{$class}{$bin}++;
            }
        }
        $total_N{$type}{$class}  = keys %{$total_N_hash{$type}{$class}};
        print "##$class $type , no filter, $total_N{$type}{$class}\n";
        $total_N_filter{$type}{$class}  = keys %{$total_N_hash_filter{$type}{$class}};
        print "##$class $type , with filter, $total_N_filter{$type}{$class}\n";
    }
}

$intersect_size = keys %intersect;

print "total number of intersected libraries: $intersect_size\n";

%size1_hash;
%size1_all_hash;

my $R = Statistics::R->new();
# $R->set('PW',$lib_cutoff);
$R->run( qq`P_mass <- function(N,m,n,x) exp(lchoose(N,x)+lchoose(N-x,m-x)+lchoose(N-m,n-x)-lchoose(N,m)-lchoose(N,n))` );
# $R->run( qq`p_val <- function(num, m,n, N){sum(sapply(num:min(m,n), P_mass, m=m, n=n, N=N))}` );
$R->run( qq`TROM <- function(N,m,n,x)  sum(sapply(x:min(m,n), P_mass, m=m, n=n, N=N)) `);

foreach my $class ("hypo","hyper")
{
    foreach my $type ("CG","CHG","CHH")
    {
        open($type, ">DOM_cutoff$cutoff\_min$min_sample_size\_$class\_$type\_DMR_overlap.txt");
        $array[0][0] = $type;
        foreach my $i (1..$n){
                foreach my $j ($i..$n){
                    if ($i == $j) {$array[$j][$i] = "NA"; next;}

                    # $size1_all = keys %{$mutant_hash{$array[$i][0]}{$type}{$class}};
                    # $size1_all_hash{$array[$i][0]} = $size1_all;
                    $size1 = keys %{$intersect{$array[$i][0]}{$type}{$class}};
                    # $size1_hash{$array[$i][0]} = $size1;

                    # $size2_all = keys %{$mutant_hash{$array[0][$j]}{$type}{$class}};
                    $size2 = keys %{$intersect{$array[0][$j]}{$type}{$class}};
                    my $overplap_count = 0;
                    foreach (keys %{$intersect{$array[$i][0]}{$type}{$class}}){
                        if (exists $intersect{$array[0][$j]}{$type}{$class}{$_}){
                            $overlap{$i}{$j}{$type}{$class}{$_} = $intersect{$array[$i][0]}{$type}{$class}{$_}; 
                        }
                        
                    }
                    $overplap_count = keys %{$overlap{$i}{$j}{$type}{$class}};

                    my $TROM_Score;
                    if ($size1 >= $min_sample_size) {
                        $R->run(qq`P_value <- -log((TROM ($total_N{$type}{$class}, $size1, $size2, $overplap_count)) + 10^(-300),base=10)`);
                        # $R->run(qq`P_value <- TROM ($total_N_filter{$type}{$class}, $size1, $size2, $overplap_count)`);
                        $TROM_Score = $R->get('P_value');
                        $array[$i][$j] = $TROM_Score;
                        # print "$type, $class, $array[$i][0], $array[0][$j], $total_N_filter{$type}{$class}, $size1, $size2, $overplap_count, $TROM_Score\n";

                    }
                    else {$array[$i][$j] = "na";}
                    if ($size2 >= $min_sample_size) {
                        $array[$j][$i] = $TROM_Score;
                    }
                    else {$array[$j][$i] = "na";}
                }
        }

        foreach my $row (@array) {
            my $line = join ("\t", @$row);
            # print @$row
            print $type "$line\n";
        }
    }
}




