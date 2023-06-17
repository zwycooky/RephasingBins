#!/usr/bin/perl

use strict;

mkdir "pbs_co_scripts";
chdir "pbs_co_scripts";

my $dir = "/public/home/hjqiu/people/zwy/Fuding_population/04phasing_FD/03phasing";

my @sca = glob "$dir/hapOut/*hap.txt";
foreach (@sca) {
my $prefix = (split /\./,(split /\//,$_)[-1])[0];
my $pbs = "#PBS -N $prefix
#PBS -l nodes=1:ppn=1

module load R/4.0.0
cd $dir/pbs_co_scripts

Rscript $dir/scripts/Hapi_co.R $_ $dir/hapOut/$prefix.co.txt

";
    open PBS,'>',"$prefix.co.pbs";
    print PBS "$pbs\n";
    close PBS;
    `qsub $prefix.co.pbs`;
    print "$prefix.co submited!\n";

}

