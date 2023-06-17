#!/usr/bin/perl

use strict;
use Cwd;
use Cwd 'abs_path';

my ($input_dir,$tmp_output,$output,$genome) = @ARGV[0,1,2,3];
my $Usage = "\n\t$0 <input bam dir> <tmp output dir> <output dir> <genome>
\n";
die $Usage unless (@ARGV == 4);

$input_dir = abs_path($input_dir);
$tmp_output = abs_path($tmp_output);
$output = abs_path($output);
$genome = abs_path($genome);
my $gatk = "gatk";
my $postfix = "sorted.marked.duplicates.bam";

## find bam file ##
chomp(my @bam = `find $input_dir -name "*.$postfix"`);
if (@bam == 0) {
	die "No bam file found:$!";
}

mkdir $tmp_output;
mkdir $output;

## get chromosome name and length ##
chomp(my @sca_id = `grep '>' $genome`);
open GENOME,'<',"$genome" or die;
my ($gname,%glen);
while (<GENOME>) {
	chomp;
	if (/>/) {
		$gname = (split />/,(split)[0])[1];
		$glen{$gname} = '';
	}else{
		$glen{$gname} += length($_);
	}
}
close GENOME;

## get bam file input command ##
my @bam_com;
foreach (@bam) {
	my $tmp = "-I $_";
	push @bam_com,$tmp;
}
my $bam_com = join(" ",@bam_com);

## Run GATK for each chromosome ##

open OUT2,'>',"Raw_split_vcf.list";
open OUT3,'>',"LSF_directory.list";

chomp(my $pwd = `pwd`);
foreach (@sca_id) {
	my $id = (split />/,(split)[0])[1];
	chdir "$pwd";
	mkdir "$tmp_output/$id";
	chdir "$tmp_output/$id";
	my $glen = $glen{$id};
	my $window = 10000000;
	my $step = 9900000;
	my $num = int($glen/$step) + 1;
	foreach (1..$num) {
		my $postfix = $_;
		my $split_s = 1 + ($_-1) * $step;
		my $split_e = $split_s + $window;
		if ($split_e >= $glen) {
			$split_e = $glen;
		}
		my $pbs_script = &lsf($id,$split_s,$split_e,$postfix);
			
		open OUT,'>',"$tmp_output/$id/$id.$postfix.pbs";
                print OUT "$pbs_script";
                close OUT;
		
		print OUT2 "$tmp_output/$id.$postfix.raw.vcf";
		print OUT3 "$tmp_output/$id/#$id.$postfix.pbs\n";

		`qsub $tmp_output/$id/$id.$postfix.pbs`;			

		if ($split_e == $glen) {
			last;
		}
	}
}

close OUT2;
close OUT3;


sub lsf {
	my ($id,$s,$e,$postfix) = @_;
	my $lsf;
	if ($postfix != 0) {
	$lsf = "#PBS -N GATK.$id.$postfix
#PBS -l nodes=1:ppn=1

module load GATK/4.3.0.0

cd $tmp_output/$id

$gatk HaplotypeCaller -R $genome $bam_com --verbosity ERROR -L $id:$s-$e -O $id.$postfix.raw.vcf
";
	}else{
	$lsf = "#PBS -N GATK.$id
#PBS -l nodes=1:ppn=1

cd $tmp_output/$id

$gatk HaplotypeCaller -R $genome $bam_com --verbosity ERROR -L $id -O $id.raw.vcf
";
	}
	return ($lsf);
}
