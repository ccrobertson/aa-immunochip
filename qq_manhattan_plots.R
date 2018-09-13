##Manhattan plots and Q-Q plots in R
#install.packages("qqman", repos="http://cran.us.r-project.org")
library("qqman") #vignette("qqman") #get help for this package


setwd("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/assoc_results")
df_raw = read.table("plink.logistic.2PC.assoc.logistic", header=TRUE)
#df_raw = read.table("plink.logistic.4PC.assoc.logistic", header=TRUE)
df_raw_gen = df_raw[df_raw$TEST=="ADD",] #remove covariate p-values
df_adjust = read.table("plink.logistic.2PC.assoc.logistic.adjusted", header=TRUE)
df_gen = merge(df_raw_gen, df_adjust, by=c("CHR", "SNP"))
dim(df_gen)



#check for inflation
qq(df_gen$P)
qq(df_gen$GC)

#exclude HLA
df_gen_no_hla = df_gen[!(df_gen$CHR==6 & df_gen$BP>29500000 & df_gen$BP<33500000),]
dim(df_gen_no_hla)
qq(df_gen_no_hla$P)

#exclude HLA & INS
df_gen_no_hla_ins = df_gen[!(df_gen$CHR==6 & df_gen$BP>29500000 & df_gen$BP<33500000) & !(df_gen$CHR==11 & df_gen$BP>1981009 & df_gen$BP<2381009),]
dim(df_gen_no_hla_ins)
qq(df_gen_no_hla_ins$P)


#manhattan plot
manhattan(df_gen, ylim = c(0,15), suggestiveline=FALSE, p="GC")
bonferroni_thresh = -log10(0.05/dim(df_gen)[1])
abline(h=bonferroni_thresh, col="blue")


#write to file 
#-- metal format
write.table(df_gen[,c("SNP", "GC")], file="results_2PCs_metalformat.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#-- tabix format?
write.table(df_gen[,c("CHR", "BP", "SNP", "OR", "GC")], file="results_2PCs_tabixformat.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#-- epacts format
epacts_format = df_gen[,c("CHR", "BP", "BP", "SNP", "GC")]
names(epacts_format) = c("#CHROM", "BEGIN", "END", "MARKER_ID", "PVALUE")
write.table(epacts_format, file="results_2PCs_epactsformat.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


######Get top hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_gen_sig = df_gen[df_gen$P<3e-7 & !is.na(df_gen$P),]
table(df_gen_sig$CHR)

#CHROM 6
chr6 = df_gen_sig[df_gen_sig$CHR==6,]
chr6[chr6$P==min(chr6$P),]

#CHROM 11
df_gen_sig[df_gen_sig$CHR==11,]
#imm_11_2141424

#CHROM 17
df_gen_sig[df_gen_sig$CHR==17,]


##Allele frequency distribution
freq = read.table("plink.frq", comment.char="~", header=TRUE)
hist(freq$MAF)
