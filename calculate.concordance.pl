#!/bin/perl/ 
use warnings;
use strict;

### THIS SCRIPT CALCULATES CONCORDANCE FOR EACH CLASSICAL HLA GENE
### The calculation is based on the formula provided in the original
### SNP2HLA paper (Jia et al 2013, Plos One) defining a formula for Acc(L) = "imputation accuracy at HLA locus L"

### INPUT: This script requires three input files
### 1. A plink format .fam file listing all the samples in your data set (see example AA_IMPUTED_validation.fam)
### 2. A dosage file as returned by SNP2HLA (one row per variant, columns = variant_id, ref_allele, alt_allele, subject1, ..., subjectN) (see example AA_IMPUTED_validation.dosage)
### 3. A ped file containing the true HLA alleles (see example T1DGC_HLA_alleles_new.ped)
### 4. A list of subjects to use for the accuracy estimate (see example AA_icaa_validation_subjects.txt) -- you may want to use all your samples, in which case this file could just be the first two columns as your .fam file described above

### OUTPUT: This script generates one output file (see example validation_alleles_dosage_new.txt) as well as prints
### an accuracy calculation for each of the 8 HLA genes to standard output.



### Get dosage dictionary (marker --> AA_icaa dosages)
my $fam = "AA_IMPUTED_validation.fam";
my @samples;
open(FAM, "<$fam") or die("Cannot read $fam\n");
while(<FAM>) {
	chomp;
	my @cols = split(/\s+/, $_);
	my $id = "$cols[0]:$cols[1]";
	push @samples, $id;
}
close(FAM);
#print scalar @samples,"\n";


my $dosage = "AA_IMPUTED_validation.dosage";
my $markerToDosage = {};
open(DOS, "grep HLA $dosage |") or die("Cannot open $dosage\n");
while(<DOS>) {
	chomp;
	my @cols = split("\t", $_);
	exit 0 if $#cols ne ($#samples+3); 
	my ($marker, $ref, $alt) = @cols[0..2]; #print join("\t", $marker, $ref, $alt),"\n";	
	for my $i(0..$#samples) {
		$markerToDosage -> {$marker} -> {$samples[$i]} = $cols[$i+3];
	}
}
close(DOS);
#print join("\t", keys %$markerToDosage),"\n";
#print join("\t", keys %{$markerToDosage->{'HLA_B_38'}}),"\n";

#Spot check for correct mapping
#print $markerToDosage->{'HLA_A_2402'}->{'46185110:46185110'},"\n";
#print $markerToDosage->{'HLA_A_101'}->{'100552:10055201'},"\n";
#print $markerToDosage->{'HLA_A_0101'}->{'CNTL5114:CNTL5114'},"\n";
#print $markerToDosage->{'HLA_A_02'}->{'9961#60053:9961#60053'},"\n";
#print $markerToDosage->{'HLA_A_02'}->{'47907511:47907511'},"\n";
#print $markerToDosage->{'HLA_A_02'}->{'CNTL4259:CNTL4259'},"\n";




### Get HLA types
my $val_subids = "AA_icaa_validation_subjects.txt";
my %validation_subjects;
open(VAL,"<$val_subids") or die("Cannot read $val_subids\n");
while(<VAL>) {
	chomp;
	my ($fid, $iid) = split(/\s+/, $_);
	$validation_subjects{$fid.":".$iid} = 1;
}
close(VAL);
#print keys %validation_subjects,"\n";

##issues with parsing ped file --> read into R and re-write to txt file
#df <- read.table("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_imputation/mkreference/T1DGC_HLA_alleles.ped", header=FALSE)
#write.table(df, "/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_imputation/mkreference/T1DGC_HLA_alleles_new.ped", quote=FALSE, row.names=FALSE, col.names=FALSE)


#initialize
my @loci = qw(HLA_A HLA_B HLA_C HLA_DPA1 HLA_DPB1 HLA_DQA1 HLA_DQB1 HLA_DRB1);
my $count = {};
my $combined_dosage = {};
for my $locus(@loci) {
	$count -> {$locus} = 0;
	$combined_dosage -> {$locus} -> {'2D'} = 0;
	$combined_dosage -> {$locus} -> {'4D'} = 0;
}

my $types = "T1DGC_HLA_alleles_new.ped";
my $outfile = "validation_alleles_dosage_new2.txt";
open(TYP, "<$types") or die("Cannot open $types\n");
open(OUT, ">$outfile") or die("Cannot write to $outfile\n");
print OUT join("\t", "FID", "IID", "ALLELE1", "ALLELE2", "DOSE1", "DOSE2", "DOSE_COMBINED"),"\n";
my $n = 0;
while(<TYP>) {
	chomp;
	my @cols = split(/\s+/, $_);
	my ($famid, $analytic_id) = @cols[0..1];
	my $id = $famid.":".$analytic_id ;	
	if (defined($validation_subjects{$id})) {	
		$n++;			
		for my $j (0..$#loci) {
			
			my $locus = $loci[$j];
			my $a1_4d = sprintf("%04d",$cols[6+2*$j]);
			my $a2_4d = sprintf("%04d",$cols[7+2*$j]); 			
			my $a1_2d = substr($a1_4d, 0, 2);
			my $a2_2d = substr($a2_4d, 0, 2);
			
			my $dosage_2d = getDoseSum($a1_2d, $a2_2d, $locus, $id);	
			my $dosage_4d = getDoseSum($a1_4d, $a2_4d, $locus, $id);
			#print join("\t", $id, $locus, $a1_2d, $a1_4d, $a2_2d, $a2_4d, $dosage_2d, $dosage_4d), "\n";
			
			$combined_dosage -> {$locus} -> {'2D'} += $dosage_2d;
			$combined_dosage -> {$locus} -> {'4D'} += $dosage_4d;
		}		
	}
}


print join("\t", "Locus", "2D_accuraccy", "4D_accuracy"),"\n";
for my $locus(@loci) {
	my $acc_2D = ($combined_dosage -> {$locus} -> {'2D'})/(2*$n);
	my $acc_4D = ($combined_dosage -> {$locus} -> {'4D'})/(2*$n);
	print join("\t", $locus, sprintf('%.3f',$acc_2D), sprintf('%.3f',$acc_4D)),"\n";
}


sub getDoseSum {
	my ($a1, $a2, $locus, $id) = @_;
	my $dose1;
	my $dose2;
	my $dosage;
	if (defined($markerToDosage->{$locus.'_'.$a1}->{$id}) && defined($markerToDosage -> {$locus.'_'.$a2} -> {$id})) {
		$dose1 = $markerToDosage->{$locus.'_'.$a1}->{$id};
		$dose2 = $markerToDosage->{$locus.'_'.$a2}->{$id};
		if ($a1 eq $a2) {
			$dosage = $dose1;
		} else {
			$dosage = $dose1 + $dose2;
		}				
	} elsif (defined($markerToDosage->{$locus.'_'.$a1}->{$id})) { #<-- in this case, a2 is not found in the reference panel, so no dosage value has been assigned to it. we set dosage for a2 = 0
		$dose1 = $markerToDosage->{$locus.'_'.$a1}->{$id};
		$dosage = $dose1;
	} elsif (defined($markerToDosage->{$locus.'_'.$a2}->{$id})) {
		$dose2 = $markerToDosage->{$locus.'_'.$a2}->{$id};
		$dosage = $dose2;
	} else {   #<-- in this case, neither allele is found in the reference panel, so no dosage value has been assigned to it. we set dosage of true alleles = 0  
		$dosage = 0;
	}
	print join("\t", $id, $locus.'_'.$a1, $locus.'_'.$a2, $dosage),"\n";
	return $dosage;
}

 










