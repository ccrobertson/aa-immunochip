#!/bin/perl
use strict;
use warnings;
use Getopt::Long;

#perl get_haplotypes.pl --in derived_data/phased_2d.txt --out derived_data/haplotypes.2d.txt --fam ../hla_imputation/imputation_results/AA_IMPUTED.fam
#perl get_haplotypes.pl --in derived_data/phased_4d.txt --out derived_data/haplotypes.4d.txt --fam ../hla_imputation/imputation_results/AA_IMPUTED.fam

my $infile;
my $outfile;
my $fam;
GetOptions("in=s" => \$infile, "out=s" => \$outfile, "fam=s" => \$fam) or die("Error in Getopt specification");

### get case control status
#my $fam = "../hla_imputation/imputation_results/AA_IMPUTED.fam";
open(FAM, "<$fam") or die("Cannot open $fam\n");
my $uniqIDtoCaseStatus = {};
while(<FAM>) {
	chomp;
	my ($FID, $IID, $FAT, $MOT, $SEX, $CASE) = split(/\s+/, $_);
	my $uniqID = join("_", $FID, $IID);
	$uniqIDtoCaseStatus -> {$uniqID} = $CASE; 
}
close(FAM);


#my $infile = "derived_data/phased_2d.txt";
open(IN, "<$infile") or die("Cannot open $infile\n");
my $colindexToFID = {};
my $header1 = <IN>;
my @cols1 = split(/\s+/, $header1);
for my $i (2..$#cols1) {
	$colindexToFID -> {$i} = $cols1[$i];
}
my $colindexToIID = {};
my $colindexToUniqID = {};
my $header2 = <IN>;
my @cols2 = split(/\s+/, $header2);
for my $i (2..$#cols2) {
	$colindexToIID -> {$i} = $cols2[$i];
	$colindexToUniqID -> {$i} = $colindexToFID -> {$i}."_".$colindexToIID -> {$i};
}
close(IN);
#print join("\t", values %$colindexToFID);
#print join("\t", values %$colindexToUniqID); 

my $uniqIDtoHaplotype = {};
open(IN, "<$infile") or die("Cannot open $infile\n");
while(<IN>) {
	chomp;
	if ($_!~/^M/) {
		next;
	}
	my @cols = split(/\s+/, $_);
	my $variant = $cols[1];
	$variant =~ /HLA_([^_]+)_.+/; my $gene = $1;	
	$variant =~ /HLA_[^_]+_([\d]+)/; my $allele = $1;
	
	for my $i (2..$#cols) {

		my $uniqID = $colindexToUniqID -> {$i};
		if ($i % 2 eq 1) {
			if ($cols[$i] eq "P") {
				$uniqIDtoHaplotype -> {$uniqID} -> {$gene} -> {'HAP1'} = $allele ;	
			}	  
		}
		if ($i % 2 eq 0) {
			if ($cols[$i] eq "P") {
				$uniqIDtoHaplotype -> {$uniqID} -> {$gene} -> {'HAP2'}  = $allele ;	
			}  
		}
	}				
}
close(IN);

#Get DRB1-DQA1-DQB1 haplotypes
#my $outfile = "derived_data/haplotypes.2d.txt"; 
open(OUT, ">$outfile") or die("Cannot open $outfile\n");
print OUT join("\t", "uniq_id", "case", "haplotype"),"\n";
#for my $uniqID (keys %$uniqIDtoCaseStatus) {
for my $uniqID (keys %$uniqIDtoHaplotype) {
	#print $uniqID,"\t";
	#print join("\t", $uniqID, $uniqIDtoCaseStatus -> {$uniqID}) ,"\n";
	
	my $hap1;
	if (defined($uniqIDtoHaplotype -> {$uniqID} -> {'DRB1'} -> {'HAP1'}) && defined($uniqIDtoHaplotype -> {$uniqID} -> {'DQA1'} -> {'HAP1'}) && defined($uniqIDtoHaplotype -> {$uniqID} -> {'DQB1'} -> {'HAP1'})) {
		$hap1 = join(":", $uniqIDtoHaplotype -> {$uniqID} -> {'DRB1'} -> {'HAP1'}, $uniqIDtoHaplotype -> {$uniqID} -> {'DQA1'} -> {'HAP1'}, $uniqIDtoHaplotype -> {$uniqID} -> {'DQB1'} -> {'HAP1'});
	} else {
		$hap1 = "NA";
	}
	
	my $hap2;
	if (defined($uniqIDtoHaplotype -> {$uniqID} -> {'DRB1'} -> {'HAP2'}) && defined($uniqIDtoHaplotype -> {$uniqID} -> {'DQA1'} -> {'HAP2'}) && defined($uniqIDtoHaplotype -> {$uniqID} -> {'DQB1'} -> {'HAP2'})) {
		$hap2 = join(":", $uniqIDtoHaplotype -> {$uniqID} -> {'DRB1'} -> {'HAP2'}, $uniqIDtoHaplotype -> {$uniqID} -> {'DQA1'} -> {'HAP2'}, $uniqIDtoHaplotype -> {$uniqID} -> {'DQB1'} -> {'HAP2'});
	} else {
		$hap2 = "NA";
	}
	
	print OUT join("\t", $uniqID, $uniqIDtoCaseStatus -> {$uniqID}, $hap1) ,"\n";
	print OUT join("\t", $uniqID, $uniqIDtoCaseStatus -> {$uniqID}, $hap2) ,"\n";
}










#X7036_X7036	1	NA
#9961#30067_9961#30067	1	NA
#X2501_X2501	1	NA
#44458711_44458711	1	NA
#44458711_44458711	1	NA
#48055811_48055811	1	NA
#CNTL4905_CNTL4905	1	NA
#04A-S-0375#0100_04A-S-0375#0100	1	NA
#7042_7042	1	NA
#46661611_46661611	1	NA
#CNTL4282_CNTL4282	1	NA
#40146811_40146811	1	NA
#2519_2519	1	NA
#40861610_40861610	2	NA
#45620611_45620611	1	NA
#48389710_48389710	2	NA
#48977711_48977711	1	NA
#CNTL4966_CNTL4966	1	NA
#44345611_44345611	1	NA
#CNTL4902_CNTL4902	1	NA
#41741011_41741011	1	NA
#21-S-0019#0100_21-S-0019#0100	1	NA
#04A-S-0386#0100_04A-S-0386#0100	1	NA
#1530_1530	1	NA
#CNTL4907_CNTL4907	1	NA
#48759211_48759211	1	NA
#CNTL38_CNTL38	1	NA
#43258911_43258911	1	NA
#CNTL5077_CNTL5077	1	NA
#9961#30122_9961#30122	1	NA
#49925711_49925711	1	NA
#9961#30048_9961#30048	1	NA
#42468510_42468510	2	NA





