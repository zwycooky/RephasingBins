#!/usr/bin/perl

use strict;

chomp(my $dir = `pwd`);
my $phasing_script = "$dir/hapi_scripts/contig_Hapi_phasing.pl";

my $Job_num = @ARGV[0];
my $Usage = "\n\t$0 <Job number>
\n";
die $Usage unless (@ARGV == 1);

mkdir "pbs_scripts";
chdir "pbs_scripts";

my @file = glob "$dir/split_vcf/*.vcf";

my ($counter,@tmp,$jobs);

foreach (@file) {
	$counter ++;
	push @tmp,$_;
	if ($counter == $Job_num) {
		my $vcf_list = join(",",@tmp);
		$jobs ++;
		$counter = 0;
		&phasing_pbs($vcf_list,$jobs);
		@tmp = ();
	}
}


sub phasing_pbs {
	my ($vcf_list,$job_id) = @_;
	my $pbs = "#PBS -N $job_id
#PBS -l nodes=1:ppn=5

#module load R/4.0.0

cd $dir
perl $phasing_script $vcf_list
";
	open PBS,'>',"$job_id.pbs";
	print PBS "$pbs";
	close PBS;
	`qsub $job_id.pbs`;
}

