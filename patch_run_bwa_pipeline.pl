#!/usr/bin/perl

## 2019-7-10 zwy ##

use strict;
use Cwd;
use Cwd 'abs_path';

my ($dir,$outdir,$genome,$bwa_index) = @ARGV[0,1,2,3];
my $Usage = "\n\t$0 <input dir of fastq data> <output dir> <reference genome> <bwa index of reference>
\n";
die $Usage unless (@ARGV == 4);

if (!-e $outdir) {
	mkdir $outdir;
}

$dir = abs_path($dir);
$outdir = abs_path($outdir);
$genome = abs_path($genome);
$bwa_index = abs_path($bwa_index);

## main scripts ##
my $scripts = "bwa-tea.pl";

## --------------------------- ##
chomp(my @input = glob "$dir/*.gz");

my $counter = 0;
my @tmp;
foreach (@input) {
	$counter ++;
	if ($counter == 2) {
		push @tmp,$_;
		my $fq1 = abs_path($tmp[0]);
		my $fq2 = abs_path($tmp[1]);
		
		my $prefix = (split /_/,(split /\//,$_)[-1])[0];
		
		chdir $outdir;
		mkdir $prefix;
		chdir $prefix;
		chomp(my $pwddir = `pwd`);
		
		## -- submit -- ##
		#print "$prefix,$pwddir,$fq1,$fq2,$genome,$bwa_index\n";
		&pbs($prefix,$pwddir,$fq1,$fq2,$genome,$bwa_index,$scripts);		

		@tmp = ();
		$counter = 0;
	}else{
		push @tmp,$_;
	}
}


sub pbs {
	my ($sample,$pwddir,$fq1,$fq2,$genome,$bwa_index,$scripts) = @_;
	my $pbsscripts = "
#BSUB -J sample\_bwa
#BSUB -n 5
#BSUB -R span[hosts=1]
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q normal

module load BWA/0.7.17
module load SAMtools/1.9

cd $pwddir

perl $scripts -1 $fq1 -2 $fq2 -f $genome -b $bwa_index -p 5 -o $sample
";
	open OUT,'>',"bwa\_$sample.lsf";
	print OUT "$pbsscripts";
	close OUT;

	`bsub <bwa\_$sample.lsf`;
	print "$sample submit!\n";

}









