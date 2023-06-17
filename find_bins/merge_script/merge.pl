#!/usr/bin/perl

use strict;
my $postfix = @ARGV[0];
my $Usage = "\n\t$0 <postfix>
\n";
die $Usage unless (@ARGV == 1);

my @file = glob "*$postfix";
my $fir_file = 0;
foreach (sort { (split /ptg000/,(split /\./,$a)[0])[1] <=> (split /ptg000/,(split /\./,$b)[0])[1] } @file) {
	if ($fir_file == 0) {
		open IN,'<',"$_" or die;
		while (<IN>) {
			chomp;
			print "$_\n";
		}
		close IN;
		$fir_file = 1;
	}else{
		my $fir = 0;
		open IN,'<',"$_" or die;
		while (<IN>) {
			chomp;
			if ($fir == 0) {
				$fir = 1;
			}else{
				print "$_\n";
			}
		}
		close IN;
	}
}

