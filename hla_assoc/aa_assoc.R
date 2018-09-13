options(stringsAsFactors=FALSE)
library(ggplot2, lib.loc="/home/ccr5ju/R/x86_64-redhat-linux-gnu-library/3.3")


#get data
setwd("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_assoc/assoc")
d = read.table("hla_assoc_aa.txt", header=TRUE)
hu = read.csv("../sources/ng.2015_Hu.supptable2.reformatted.csv")

#parse MEGA data
snp_split = strsplit(d$SNP, split="_")
d$Gene = sapply(snp_split,"[[", 2)
d$aa_pos = sapply(snp_split,"[[", 3)
d$aa_pos_numeric=d$aa_pos
d$aa_pos_numeric[grep("x", d$aa_pos)] = sapply(strsplit(d$aa_pos[grep("x", d$aa_pos)],split="x"),"[[",1)
d$aa_pos_numeric = as.numeric(d$aa_pos_numeric)
d$bp_pos = as.numeric(sapply(snp_split,"[[", 4))
#d$allele = sapply(snp_split,"[[", 5)
d$Ancestry="African"
d$AF = d$FRQ_A
d$AF_AA = d$AF
d$empOR = (d$FRQ_A/(1-d$FRQ_A))/(d$FRQ_U/(1-d$FRQ_U))

#parse Hu data
hu$SNP = hu$SNP2HLA.variant.name
hu$OR = as.numeric(hu$OR.Total.)
hu$OR_EUR = hu$OR
hu$P = hu$log10P.Total.
hu$P_EUR = hu$P
hu$aa_pos = hu$Amino.acid.position
hu_aa_1 = hu[grep("AA_", hu$SNP2HLA.variant.name),]
hu_aa = hu_aa_1[hu_aa_1$Case.frequency>0.001 | hu_aa_1$Control.frequency>0.001,]
hu_aa$aa_pos_numeric=as.numeric(hu_aa$aa_pos)
hu_aa$Ancestry="European"
hu_aa$AF = hu_aa$Case.frequency
hu_aa$AF_EUR = hu_aa$AF

#merge datasets
m_short = merge(d[,c("Gene","SNP","aa_pos_numeric","OR", "AF_AA", "P")], hu_aa[,c("Gene","SNP","aa_pos_numeric","OR_EUR","AF_EUR", "P_EUR")], by=c("Gene","SNP", "aa_pos_numeric"))
m_short$ratio = m_short$OR/m_short$OR_EUR
m_long = rbind(d[,c("Gene","SNP","aa_pos_numeric","OR", "AF","Ancestry")], hu_aa[,c("Gene","SNP", "aa_pos_numeric","OR","AF","Ancestry")])
m_long$jointEffect = m_long$OR*m_long$AF

plot(-log10(m_short$P),-m_short$P_EUR)

#ggplot(m_long[m_long$Gene=="A",], aes(x=aa_pos_numeric, y=log(OR))) + geom_point(aes(color=Ancestry))
#ggplot(m_long[m_long$Gene=="C",], aes(x=aa_pos_numeric, y=log(OR))) + geom_point(aes(color=Ancestry)) + geom_hline(aes(yintercept=0)) + ggtitle("C") + xlab("Amino Acid Residue") + facet_grid(Ancestry~.)
#theme(aspect.ratio=4/3)
getPlots = function(gene) {
	axisstyle <- element_text(face = "bold", color = "blue", size = 18)
	p1 = ggplot(m_long[m_long$Gene==gene,], aes(x=aa_pos_numeric, y=log(OR))) + geom_point(aes(color=Ancestry)) + 
			geom_hline(aes(yintercept=0)) + 
			ggtitle(gene) + 
			#xlab("Amino Acid Residue") + 
			xlab("") + ylab("") +
			facet_grid(Ancestry~.) +
			theme(axis.text = axisstyle) +
			theme(aspect.ratio=1/6)
	print(p1)
	p2 = ggplot(m_long[m_long$Gene==gene,], aes(x=aa_pos_numeric, y=AF)) + geom_point(aes(color=Ancestry)) + 
			geom_hline(aes(yintercept=0)) + 
			ggtitle(gene) + 
			#xlab("Amino Acid Residue") + 
			xlab("") + ylab("") +
			facet_grid(Ancestry~.) +
			theme(axis.text = axisstyle) +
			theme(aspect.ratio=1/6)
	print(p2)
	#p2 = ggplot(m_short[m_short$Gene==gene,], aes(x=aa_pos_numeric, y=log(ratio))) + geom_point() + geom_hline(aes(yintercept=0)) + ggtitle(gene) + xlab("Amino Acid Residue")
	#print(p2)
	#p3 = ggplot(m_long[m_long$Gene==gene,], aes(x=AF, y=log(OR))) + geom_point(aes(colour=Ancestry)) + ggtitle(gene)
	#print(p3)
	#p4 = ggplot(m_long[m_long$Gene==gene,], aes(x=aa_pos_numeric, y=log(jointEffect))) + geom_point(aes(colour=Ancestry)) + ggtitle(gene)
	#print(p4)
}

pdf("aa_association_eur_vs_aa.pdf")
genes = c("A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1")
for ( i in 1:length(genes)) {	
	getPlots(genes[i])
}
dev.off()
#scp aa_association_eur_vs_aa.pdf "$local/AA_Immunochip"




###Looking at specific sites
b1 = m_short[m_short$Gene=="B",c("SNP","OR","OR_EUR")]
b1_signal=b1[b1$OR_EUR<0.15,]
d[d$SNP%in% b1_signal$SNP,]