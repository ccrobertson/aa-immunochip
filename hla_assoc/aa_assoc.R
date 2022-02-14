options(stringsAsFactors=FALSE)
library(ggplot2)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Get aminoa acid summary stats from AFR and EUR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/data1/ccr5ju/AA_Immuno/hla_assoc/assoc")

### get icaa (AFR) summary stats
#d = read.table("hla_assoc_aa.txt", header=TRUE)
load(file="conditional_analysis_object.aa.RData")
d = data.frame(out.aa[[1]][[5]])
d$SNP = row.names(d)
d$FRQ_A = d$case_AF
d$FRQ_U = d$control_AF

# parse icaa data
snp_split = strsplit(d$SNP, split="_")
d$Gene = sapply(snp_split,"[[", 2)
d$aa_pos = sapply(snp_split,"[[", 3)
d$aa_pos_numeric=d$aa_pos
d$aa_pos_numeric[grep("x", d$aa_pos)] = sapply(strsplit(d$aa_pos[grep("x", d$aa_pos)],split="x"),"[[",1)
d$aa_pos_numeric = as.numeric(d$aa_pos_numeric)
d$bp_pos = as.numeric(sapply(snp_split,"[[", 4))
#d$allele = sapply(snp_split,"[[", 5)
d$Ancestry="African"
d$case_AF_AA = d$FRQ_A
d$unadjustedOR = (d$FRQ_A/(1-d$FRQ_A))/(d$FRQ_U/(1-d$FRQ_U))
#d$empOR = (d$FRQ_A/(1-d$FRQ_A))/(d$FRQ_U/(1-d$FRQ_U))
#d$unadjustedOR = (d$case_AF/(1-d$case_AF))/(d$control_AF/(1-d$control_AF))

# filter by MAF
d$case_MAF = sapply(d$case_AF, function(x) {min(x, 1-x)})
d$control_MAF = sapply(d$control_AF, function(x) {min(x, 1-x)})
d = d[d$case_MAF>0.025 | d$control_MAF>0.025,]


### get EUR summary stats
hu = read.csv("../sources/ng.2015_Hu.supptable2.reformatted.csv")

# parse Hu data
hu$SNP = hu$SNP2HLA.variant.name
hu$OR = as.numeric(hu$OR.Eur.)
hu$OR_EUR = hu$OR
hu$P = hu$log10P.Eur.
hu$P_EUR = hu$P
hu$aa_pos = hu$Amino.acid.position
hu_aa_1 = hu[grep("AA_", hu$SNP2HLA.variant.name),]
hu_aa = hu_aa_1[hu_aa_1$Case.frequency>0.001 | hu_aa_1$Control.frequency>0.001,]
hu_aa$aa_pos_numeric=as.numeric(hu_aa$aa_pos)
hu_aa$Ancestry="European"
hu_aa$case_AF = hu_aa$Case.frequency
hu_aa$case_AF_EUR = hu_aa$case_AF

# filter by MAF
hu_aa$case_MAF = sapply(hu_aa$Case.frequency, function(x) {min(x, 1-x)})
hu_aa$control_MAF = sapply(hu_aa$Control.frequency, function(x) {min(x, 1-x)})
hu_aa = hu_aa[hu_aa$case_MAF>0.025 | hu_aa$control_MAF>0.025,]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge EUR and AFR summary stats
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Merge - keeping all variants from each data set
m_short = merge(d[,c("Gene","SNP","aa_pos_numeric","OR", "unadjustedOR", "case_AF_AA", "P")], hu_aa[,c("Gene","SNP","aa_pos_numeric","OR_EUR","case_AF_EUR", "P_EUR")], by=c("Gene","SNP", "aa_pos_numeric"), all=TRUE)


### Create long format
m_long = rbind(d[,c("Gene","SNP","aa_pos_numeric","OR", "case_AF","Ancestry")], hu_aa[,c("Gene","SNP", "aa_pos_numeric","OR","case_AF","Ancestry")])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot correlation between effect sizes between ancestry groups
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdf("aa_eur_vs_afr_logOR.pdf")
ggplot(m_short, aes(x=log(OR_EUR), y=log(OR))) + geom_point() +
	xlab("European ancestry log(OR)") + ylab("African ancestry log(OR)") +
	facet_wrap(~Gene, ncol=2) +
	theme_bw()

ggplot(m_short, aes(x=log(OR_EUR), y=log(unadjustedOR))) + geom_point() +
		xlab("European ancestry log(OR)") + ylab("African ancestry log(OR) - unadjusted") +
		facet_wrap(~Gene, ncol=2) +
		theme_bw()

dev.off()


getPlots = function(gene) {
	axisstyle <- element_text(face = "bold", color = "blue", size = 10)
	p1 = ggplot(m_long[m_long$Gene==gene,], aes(x=aa_pos_numeric, y=log(OR))) + geom_point(aes(color=Ancestry)) +
			geom_hline(aes(yintercept=0)) +
			ggtitle(gene) +
			#xlab("Amino Acid Residue") +
			xlab("") + ylab("log(OR)") +
			facet_grid(Ancestry~.) +
			theme_bw() +
			theme(axis.text = axisstyle) +
			theme(aspect.ratio=1/4)
	print(p1)
	p2 = ggplot(m_long[m_long$Gene==gene,], aes(x=aa_pos_numeric, y=case_AF)) + geom_point(aes(color=Ancestry)) +
			geom_hline(aes(yintercept=0)) +
			ggtitle(gene) +
			#xlab("Amino Acid Residue") +
			xlab("") + ylab("Case Allele Frequency") +
			facet_grid(Ancestry~.) +
			theme_bw() +
			theme(axis.text = element_text(face = "bold", color = "blue", size = 10)) +
			theme(aspect.ratio=1/4)
	print(p2)
}

pdf("aa_association_eur_vs_aa.pdf")
genes = c("A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1")
for ( i in 1:length(genes)) {
	getPlots(genes[i])
}
dev.off()
