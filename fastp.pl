#!/usr/bin/perl

use strict;
use Cwd 'abs_path';

my ($inputdir,$cpu) = @ARGV[0,1];
my $Usage = "\n\t$0 <input fastq dir> <threads [10]>
\n";
die $Usage unless (@ARGV == 2);

chomp(my $pwd = `pwd`);

$inputdir = abs_path($inputdir);

my $reads_dir = "clean_reads";
if (!-e $reads_dir) {
	mkdir $reads_dir;
	$reads_dir = abs_path($reads_dir);
}

my $fastp_pbs_dir = "fastp_pbs";
if (!-e $fastp_pbs_dir) {
	mkdir $fastp_pbs_dir;
}

chomp(my @fq = `find $inputdir -name "*.gz"`);
my $count = 0;
my ($fq1,$fq2);
foreach (sort @fq) {
	$count ++;
	if ($count < 2) {
		$fq1 = abs_path($_);
	}else{
		$fq2 = abs_path($_);
		my $id = (split /\//,$fq2)[-3];
		$count = 0;
		&pbs($fq1,$fq2,$id,$reads_dir,$fastp_pbs_dir);
	}
}

sub pbs {
	my ($fq1,$fq2,$id,$reads_dir,$pbs_dir) = @_;
	my $pbs = "#PBS -N $id.fastp
#PBS -l nodes=1:ppn=10

cd $reads_dir
fastp -i $fq1 -I $fq2 -o $id.R1.clean.fq.gz -O $id.R2.clean.fq.gz --thread=$cpu

";
	chdir $pwd;
	open PBS,'>',"$pbs_dir/$id.pbs";
	print PBS "$pbs\n";
	close PBS;
	chdir $pbs_dir;
	`qsub $id.pbs`;
}

