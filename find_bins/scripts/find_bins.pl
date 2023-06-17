#!/usr/bin/perl

use strict;
use Cwd 'abs_path';

my $input_dir = @ARGV[0];
my $Usage ="\n\t $0 <hapOut>
\n";
die $Usage unless (@ARGV == 1);

my @co_file = glob "$input_dir*co.txt";
my @hap_file = glob "$input_dir*hap.txt";

my %co_file;
foreach (@co_file) {
	my $file_path = abs_path($_);
	my $contig_id = (split /\./, (split /\//,$file_path)[-1])[0];
	$co_file{$contig_id} = $file_path;
}

my %hap_file;
foreach (@hap_file) {
        my $file_path = abs_path($_);
        my $contig_id = (split /\./, (split /\//,$file_path)[-1])[0];
        $hap_file{$contig_id} = $file_path;
}

my $find_bin_script = "/public/home/hjqiu/people/zwy/Fuding_population/04phasing_FD/03phasing/find_bins/scripts/find_bin.R";
my $find_bin_script2 = "/public/home/hjqiu/people/zwy/Fuding_population/04phasing_FD/03phasing/find_bins/scripts/get_co_from_sca.R";

if (!-e "results") {
	mkdir "results";
	chdir "results";
}else{
	chdir "results";
}

foreach (sort keys %hap_file) {
	my $contig_id = $_;
	if (exists $co_file{$_}) {
		!system "Rscript $find_bin_script $co_file{$_} $hap_file{$_}" or die "Error with find bins:$!";
		#print "Rscript $find_bin_script $co_file{$_} $hap_file{$_}\n";
	}else{
		!system "Rscript $find_bin_script2 $hap_file{$_} /public/home/hjqiu/people/zwy/Fuding_population/04phasing_FD/03phasing/find_bins/FD_hifi.asm.bp.p_ctg.length.txt" or die "Error with find bins2:$!";
	}
	#exit;
}


