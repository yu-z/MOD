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

$folder100bp = "/u/home/j/jxzhai/MCDB_folder/Projects/BigdataBSseq/2016/503Lib_100bp/";


# obtain the list of WT libraries and save as hash
open (WT,$wtlist);
my %WT_hash;
# my %mutant_hash;
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

open (REF,"/u/home/j/jxzhai/MCDB_folder/Projects/BigdataBSseq/2016/Lib_list/ref.lst");
my %REF_hash;

while (<REF>){
    chomp;
    next if (length($_)<1);
    my ($gsmid, $sample) = split/\t/;
    $REF_hash{$gsmid} = $sample;
}


$WT_size = keys %WT_hash;
$MU_size = keys %MU_hash;

# $cutoff = $WT_size;
# $cutoff = $WT_size * 0.5;
$cutoff = 33;

$min_sample_size = 0;



print "total number of WT libraries: $WT_size , cutoff set at $cutoff , total number of MU libraries: $MU_size , \n";

# die;

my %mutant_hash;
my %goodbin_hash;
my $goodbin_file;

opendir (DIR100, "$folder100bp") or die "Couldn't open $seqdir: $!\n";
while (defined(my $seqfile = readdir(DIR100))) {

    next if (($seqfile =~ m/^\.$/) or ($seqfile =~ m/^\.\.$/));
    next if !($seqfile =~ m/$run_type\.100\.gz$/);
    # print "$seqfile\n";

    $seqfile =~ m/(\S+)\.(\S+)\.100\.gz/;
    $name = $1;
    my $type = $2;
    $name =~ s/#\d+//;
    next if !(exists($MU_hash{$name}));
    print "##$name##\n";

    print "$seqfile\n";

    print "$name\t$type\n";

    open(FH, "gunzip -c $folder100bp$seqfile |") || die "can't open pipe to $file";
    while (<FH>) {
        chomp;
        next if (m/chr/i);
        @data = split/\t/;
        my $bin_key = "\"$data[0]\"\_$data[1]";
        # print "$bin_key\n";
        ###### count of bin_key occurrence
        if ($data[5]>=4) {
            $goodbin_hash{$name}{$type}{$bin_key}++;
        }
    }
    close FH;
    my $goodbin_number = keys %{$goodbin_hash{$name}{$type}};
    $goodbin_file++;
    print "$name\t$type\t$goodbin_number\n";
 
}
print "loaded goodbin files: $goodbin_file\n";
closedir DIR100;


opendir (DIR, $seqdir) or die "Couldn't open $seqdir: $!\n";

while (defined(my $seqfile = readdir(DIR))) {
    next if (($seqfile =~ m/^\.$/) or ($seqfile =~ m/^\.\.$/));
    next if !($seqfile =~ m/DMR/);
    next if !($seqfile =~ m/^$run_type\./);
    
    $seqfile =~ /(.*)\.DMR\.(.*)_VS_(.*)\.txt/;
    my $type = $1;
    my $lib1 = $2; $lib1 =~ s/#\d+//;
    my $lib2 = $3; $lib2 =~ s/#\d+//;
    
    if (exists $WT_hash{$lib1}) {
        # if (exists $WT_hash{$lib2}) {next;}
        next if !(exists $MU_hash{$lib2});
        $WT = $lib1; $mutant = $lib2;
        # print "$WT\t$mutant\n" if ($mutant =~ m/CG|CHG|CHH/i);
    }
    elsif (exists $WT_hash{$lib2}) {
        # if (exists $WT_hash{$lib1}) {next;}
        next if !(exists $MU_hash{$lib1});
        $WT = $lib2; $mutant = $lib1;
    }
    else {
        next;
    }
    print "$WT\t$mutant\n";

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

foreach my $lib1 (sort keys %mutant_hash){
    $n++;
    $array[$n][0] = $lib1;
    $array[0][$n] = $lib1;
}
print "total number of mutant libraries: $n\n";


my %intersect;
my %overlap;

foreach my $class ("hypo","hyper")
{
    #foreach my $type ("CG","CHG","CHH")
    my $type = $run_type;
    {
        foreach my $lib1 (sort keys %mutant_hash){
            foreach my $bin (keys %{$mutant_hash{$lib1}{$type}{$class}}){
                if ($mutant_hash{$lib1}{$type}{$class}{$bin} >= $cutoff) {
                $intersect{$lib1}{$type}{$class}{$bin} = $mutant_hash{$lib1}{$type}{$class}{$bin};
                }
            }
        }
    }
}

$intersect_size = keys %intersect;
print "total number of intersected libraries: $intersect_size\n";

%size1_hash;
%size1_all_hash;

foreach my $class ("hypo","hyper")
{
    #foreach my $type ("CG","CHG","CHH")
    my $type = $run_type;
    {
        open($type, ">cutoff$cutoff\_min$min_sample_size\_$class\_$type\_DMR_overlap.txt");
        $array[0][0] = $type;
        foreach my $i (1..$n){
                foreach my $j ($i..$n){

                    my $size1_all = keys %{$mutant_hash{$array[$i][0]}{$type}{$class}};
                    $size1_all_hash{$array[$i][0]} = $size1_all;
                    my $size1 = keys %{$intersect{$array[$i][0]}{$type}{$class}};
                    $size1_hash{$array[$i][0]} = $size1;

                    if ($i == $j) {$array[$j][$i] = "NA"; next;}

                    # $size2_all = keys %{$mutant_hash{$array[0][$j]}{$type}{$class}};
                    # $size2 = keys %{$intersect{$array[0][$j]}{$type}{$class}};
                    
                    # $array[$i][$j] = $hash{$array[$i][0]}{$type}{$class}{$bin};
                    # #column
                    # $array[$j][$i] = $hash{$type}{$array[$i][0]}{$array[0][$j]}[2];
                    # print "$i\t$j\t$array[$i][0]\t$array[0][$j]\t$array[$i][$j]\t$array[$j][$i]\n";
                    #print "$array[$j][$i]\n";
                    my $overplap_count = 0;
                    my $size1_adj = 0;
                    my $size2_adj = 0;
                    ### count size1
                    # print "$array[0][$j]\t$type\n";
                    foreach (keys %{$intersect{$array[$i][0]}{$type}{$class}}){
                        # print "###$_##\n";
                        if (exists $goodbin_hash{$array[0][$j]}{$type}{$_}){
                            # print "$_\n";
                            $size1_adj++;

                        }
                        if (exists $intersect{$array[0][$j]}{$type}{$class}{$_}){
                            $overlap{$i}{$j}{$type}{$class}{$_} = $intersect{$array[$i][0]}{$type}{$class}{$_}; 
                        }
                    }

                    ### count size2

                    foreach (keys %{$intersect{$array[0][$j]}{$type}{$class}}){
                        if (exists $goodbin_hash{$array[$i][0]}{$type}{$_}){
                            $size2_adj++;
                        }
                    }

                    $overplap_count = keys %{$overlap{$i}{$j}{$type}{$class}};
                    # my $ratio1;
                    # my $ratio2;
                    
                    if ($size1_adj >= $min_sample_size) {
                        # $ratio1 = sprintf("%.4f", $overplap_count/$size1); 
                        $array[$i][$j] = "\=$overplap_count\/$size1_adj";
                    }
                    else {$array[$i][$j] = "na";}
                    if ($size2_adj >= $min_sample_size) {
                        # $ratio2 = sprintf("%.4f", $overplap_count/$size2);
                        $array[$j][$i] = "\=$overplap_count\/$size2_adj";
                    }
                    else {$array[$j][$i] = "na";}
                    
                    # $array[$i][$j] = $ratio1;
                    # $array[$j][$i] = $ratio2;

                     # print "$type\n$array[$i][0]\t$size1_all\t$size1\n";
                     # print "$array[0][$j]\t$size2_all\t$size2\n";
                     # print "$overplap_count\t$ratio1\t$ratio2\n";
                }
            }
        #print "$array[0][0]\t$array[0][1]\t$array[0][2]\n";
        #print "$array[1][0]\t$array[1][1]\t$array[1][2]\n";
        #print "$array[2][0]\t$array[2][1]\t$array[2][2]\n";
        foreach my $row (@array) {
            my $line = join ("\t", @$row);
            $line =~ s/lib_//gi;
            # print @$row
            print $type "$size1_all_hash{$$row[0]}\t$size1_hash{$$row[0]}\t$REF_hash{$$row[0]}\t$line\n";
            # print  "$size1_all_hash{$$row[0]}\t$size1_hash{$$row[0]}\t$REF_hash{$$row[0]}\t$line\n";
            # print $type "$line\n";
        }
    }
}

