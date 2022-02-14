#!/bin/perl
use strict;
use warnings;
use Getopt::Long;


#perl extractCols.pl --file derived_data/imputed_add_geno.raw --cols FID,IID,PAT,MAT,SEX,PHENOTYPE,HLA_DRB1_1601_P,HLA_DRB1_1602_P

my $file;
my $cols_to_keep;
my $out;
GetOptions("file=s" => \$file, "cols=s" => \$cols_to_keep, "out=s" => \$out) 
or die("Error in command line arguments\nUsage: perl extractCols.pl --file <filename> --cols <colname1,colname2,colname3,... --out <outfilename>\n");


#This script is for retrieving columns from a file based on the header name
#It assumes input file is white-space delimited and contains a single header line
#It will parse the header and determine which fields to print in the rest of the file
#Produces tab-delimited file

#Parse colnames
my @keep_colnames = split(",", $cols_to_keep);
my %keepHash;
for my $colname (@keep_colnames) {
	$keepHash{$colname} = 1;
} 


#Get header
open(IN, "<$file") or die("Cannot open $file\n");
my $header = <IN>;
my @header_cols = split(/\s+/, $header);
my @keep_colnums;
foreach my $i (0..$#header_cols) {
	my $colname = $header_cols[$i];
	if (defined($keepHash{$colname})) {
		push @keep_colnums, $i
	}
}
close(IN);


#Get columns
open(IN, "<$file") or die("Cannot open $file\n");
while(<IN>) {
	chomp;
	my @cols = split(/\s+/, $_);
	print join("\t", @cols[@keep_colnums]),"\n";
}
close(IN); 
			