#!/usr/bin/perl

use strict;

chomp(my $dir = `pwd`);
my $vcf_list = @ARGV[0];
my $Usage = "\n\t$0 <vcf list>
\n";
die $Usage unless (@ARGV == 1);

my @file = (split /,/,$vcf_list);

#`module load R/4.0.0`;

foreach (@file) {
	my $contig = (split /\./,(split /\//,$_)[-1])[0];
	!system "Rscript $dir/hapi_scripts/Hapi_step_by_step_parallel.2.R $_ $dir/hapOut/$contig.hap.txt $dir" or warn "warnning with hapi phasing $contig:$!";
	#print "Rscript $dir/scripts/Hapi_step_by_step_parallel.2.R $_ $dir/hapOut/$contig.hap.txt\n";
	print "$contig DONE!\n";
}


