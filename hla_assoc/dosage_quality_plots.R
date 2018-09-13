library(ggplot2)
library(rtf)
setwd("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_assoc")


d = read.table("derived_data/dosages.4d", header=FALSE, sep="\t")
geno = d[,c(-1,-2,-3)]
d$g = sapply(strsplit(d[,1], split="_"), "[[", 2)

quantiles = list()
genes = c("A","B","C","DQA1","DQB1","DRB1","DPA1","DPB1")
df = data.frame()
for (i in 1:length(genes)) {
	gene = genes[i]
	dose_sum = colSums(geno[d$g==gene,])
	#quantiles[[i]][1] = gene
	#quantiles[[i]][2] = quantile(dose_sum, probs=c(0.025, 0.975))
	quantiles[[i]] = quantile(dose_sum, probs=c(0.025, 0.5, 0.975))
	df = rbind(df, data.frame(dose_sum, gene))
}
out = data.frame(gene=c(genes,"combined"), round(rbind(t(do.call(cbind, quantiles)), quantile(df$dose_sum, probs=c(0.025, 0.5, 0.975))), digits=2))

rtffilename = "tables/dosage_sum_quantiles.doc"
rtffile <- RTF(rtffilename)
addParagraph(rtffile, "Supplementary Table #. Quartiles for the sum of dosages across alleles by locus")
addTable(rtffile, out)
done(rtffile)

pdf("tables/dosage_quality_boxplots.pdf")
ggplot(df, aes(x=gene, y=dose_sum)) + geom_boxplot() + xlab("Locus") + ylab("Sum of Dosages")+ geom_hline(aes(yintercept=2, col="red"))
ggplot(df, aes(x=gene, y=log2(dose_sum))) + geom_boxplot() + xlab("Locus") + ylab("log2 Sum of Dosages") + geom_hline(aes(yintercept=1, col="red"))
dev.off()


#scp tables/dosage_* "$local/AA_Immunochip/tables_and_figures/"