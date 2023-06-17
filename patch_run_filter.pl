#!/usr/bin/perl

use strict;
use Cwd 'abs_path';

my ($raw_vcf_dir,$output) = @ARGV[0,1];
my $Usage = "\n\t$0 <raw vcf dir> <output>
\n";
die $Usage unless (@ARGV == 2);

$raw_vcf_dir = abs_path($raw_vcf_dir);
$output = abs_path($output);
my $gatk_filter = "/public/home/hjqiu/people/zwy/single_cell_2021/BGI_version/01phasing/01SNP_calling/final_output/scripts/GATK_hard_filter.pl";
my $get_pass = "/public/home/hjqiu/people/zwy/single_cell_2021/BGI_version/01phasing/01SNP_calling/final_output/scripts/get_PASS_SNPs.pl";

chomp(my @raw_vcf = `find $raw_vcf_dir -name "*.raw.vcf"`);

if (! -e "filter_pbs_scripts") {
	mkdir "filter_pbs_scripts";
}
my $pbs_dir = abs_path("filter_pbs_scripts");
chdir $pbs_dir;

my $job = 0;
foreach (@raw_vcf) {
	$job ++;
	my $vcf_file = abs_path($_);
	my $file_name = (split /\.raw\.vcf/,$vcf_file)[0];
	my $out_name = $file_name . ".filtered.vcf";
	my $out_PASS = $file_name . ".filtered.PASS.vcf";
	
	open PBS,'>',"$job.pbs";
	my $pbs = "#PBS -N patch
#PBS -l nodes=1:ppn=1

module load GATK/4.3.0.0
cd $pbs_dir

$gatk_filter $vcf_file $out_name
$get_pass $out_name > $out_PASS
";
	print PBS $pbs;
	close PBS;
	`qsub $job.pbs`;
}



