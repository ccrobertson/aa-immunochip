options(stringsAsFactors=FALSE)
setwd("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_assoc")
hla_2 = read.table("assoc_results/hla_assoc_2_digit.txt", header=TRUE)		
hla_4 = read.table("assoc_results/hla_assoc_4_digit.txt", header=TRUE)
aa = read.table("assoc_results/hla_assoc_aa.txt", header=TRUE)
snp = read.table("assoc_results/hla_assoc_snps.txt", header=TRUE)


getTopVar = function(df) {
	common = df[!is.na(df$P),]
	common[common$P==min(common$P),]
}
getTopVar(hla_2)
getTopVar(hla_4)
getTopVar(snp)
getTopVar(aa)


#install.packages("ggplot2", repos="https://cloud.r-project.org")
library(ggplot2)


#note: amino acid (and snp) positions are provided in ../SNP2HLA_package_v1.0.3/MakeReference/HLA_DICTIONARY_AA.map and HLA_DICTIONARY_SNPS.map
#extraced from EMBL-EBI Immunogenetics HLA Database (http://www.ebi.ac.uk/ipd/imgt/hla/)
#also see http://hla.alleles.org
#awk '$2~/^AA_/' ../hla_imputation/SNP2HLA_package_v1.0.3/MakeReference/HLA_DICTIONARY_AA.map > hla_aa.map

#get protein sequence boundaries
map = read.table("hla_aa.map")
map$gene = sapply(sapply(X=map$V2, FUN=strsplit, split="_", fixed=TRUE), function(l) l[[2]])
genes = unique(map$gene)
gene_info = list()
for (i in 1:length(genes)) {
	gene = genes[i]
	min = min(map$V4[map$gene==gene])-1e3 #widen range so that ggplot will print lines
	max = max(map$V4[map$gene==gene])+1e3
	gene_info[[i]] = list(gene, min, max)
}
gene_info


#plot association results
makePlot = function(df) {
	df$Mb = df$BP/1e6
	p = ggplot(df) + geom_point(aes(x=Mb, y=-log10(P)))  
	label_loc = 100
	add_lines <- function(i){
		annotate("pointrange", x=gene_info[[i]][[2]]/1e6, y=label_loc+5*(i-2), ymin=0, ymax=label_loc+5*(i-2), color="red")
	}
	add_labels <- function(i){
		annotate("text", x=gene_info[[i]][[2]]/1e6, y=label_loc+5*i, label=gene_info[[i]][[1]], color="red")
	}
	p + lapply(c(1:length(gene_info)), FUN=add_lines) + lapply(c(1:length(gene_info)), FUN=add_labels) 
}
makePlot(hla_2)
makePlot(hla_4)
makePlot(aa)
makePlot(snp)
