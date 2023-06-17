#!/usr/bin/perl

use strict;
use Getopt::Std;
our($opt_1,$opt_2,$opt_f,$opt_b,$opt_p,$opt_o);
getopts('1:2:f:b:p:o:');

my $left_reads = $opt_1;
my $right_reads = $opt_2;
my $reference = $opt_f;
my $bwa_index = $opt_b;
my $threads = (defined $opt_p)?$opt_p:6;
my $output = $opt_o;

my $Usage = "\n$0 -1 -2 -f <options>
	-1 <FILE>  left reads
	-2 <FILE>  right reads
	-f <FILE>  reference genome
	-b <INDEX> BWA index
	-p <INT>   number of threads [6]
	-o <STR>   output prefix
	\n";
die $Usage unless ($opt_1 && $opt_2 && $opt_f && $opt_b && $opt_o);
#######################################################
# my software #
#my $bwa = "/public/home/hyyu/zwy/software/bwa/bwa";
#my $java = "/public/home/hyyu/zwy/software/java/jre1.8.0_281/bin/java";
#my $picard = "/public/home/hyyu/zwy/software/picard.jar";
#my $samtools = "~/test/genome/samtools-1.10/samtools";
#my $gatk = "~/zwy/software/my_bin/gatk";
#my $Unique_maped = "/home01/haiji/caculation/single_cell/zwy/scripts/Unique_map_finder.pl";
my $bwa = "bwa";
my $java = "java";
my $picard = "~/software/picard.jar";
my $samtools = "samtools";


# checking the index of reference genome #
if (!-e "$reference.fai") {
	!system "$samtools faidx $reference" or die;
}
my $prefix = (split /\.fa/,$reference)[0];
if (!-e "$prefix.dict") {
	!system "$java -jar $picard CreateSequenceDictionary R=$reference O=$prefix.dict" or die;
}


# get read group #
chomp(my $header = `zcat $left_reads |head -1`);
my ($id1,$id2) = (split /:/,$header)[2,3];
my $id = "$id1.$id2";

my $sample_id = (split /_/,(split /\//,$output)[-1])[0];
!system "$bwa mem -t $threads -M -R \"\@RG\\tID:bwa\\tPL:Illumina\\tLB:1\\tSM:$sample_id\\tPU:bwa\" $bwa_index $left_reads $right_reads > $output.sam" or die "Error with bwa mem:$!";
#print "$bwa mem -t $threads -M -R \"\@RG\tID:bwa\tPL:Illumina\tLB:1\tSM:$sample_id\tPU:bwa\" $bwa_index $left_reads $right_reads > $output.sam";

## filter for unique reads ##
#!system "perl $Unique_maped $output.sam 0.8 $output.uni.sam";

## picard ##
# save the raw bam file #
!system "$java -jar $picard SortSam I=$output.sam O=$output.sorted.bam SORT_ORDER=coordinate" or die "Error with SortSam:$!";
unlink "$output.sam";

!system "$java -jar $picard MarkDuplicates I=$output.sorted.bam O=$output.sorted.marked.duplicates.bam M=$output.marked_dup_metrics.txt REMOVE_DUPLICATES=true" or die "Error with MarkDuplicates:$!";
unlink "$output.sorted.bam";

# get the unique mapped bam file for further variants calling #
#!system "$java -jar $picard SortSam I=$output.uni.sam O=$output.uni.sorted.bam SORT_ORDER=coordinate" or die "Error with SortSam:$!";
#unlink "$output.uni.sam";
#!system "$java -jar $picard MarkDuplicates I=$output.uni.sorted.bam O=$output.uni.sorted.marked.duplicates.bam M=$output.uni.marked_dup_metrics.txt REMOVE_DUPLICATES=true" or die "Error with MarkDuplicates:$!";
#unlink "$output.uni.sorted.bam";

## index bam ##
!system "$samtools index $output.sorted.marked.duplicates.bam" or die "ERROR with index:$!";
#!system "$samtools index $output.sorted.bam" or die "ERROR with index:$!";

## GATK SNP calling ##
#!system "$gatk HaplotypeCaller -R $reference -I $output.uni.sorted.marked.duplicates.bam --emit-ref-confidence GVCF -O $output.raw.gvcf.vcf" or die "ERROR with GATK HaplotypeCaller";





