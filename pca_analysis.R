#install.packages("ggplot2", repos="http://cran.us.r-project.org")
#source("http://bioconductor.org/biocLite.R")
#biocLite("CGEN") #Bhattacharjee AJHG 2010
#install.packages("scatterplot3d", repos="http://cran.us.r-project.org")

options(stringsAsFactors=FALSE)
library(ggplot2)
source("/home/ccr5ju/multiplot.R")
library("CGEN")
library(scatterplot3d)


##Read in PC and subject info files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pcs = read.table("/m/jdrfdn_scratch/users/projects/IMCHIP/20160316_African_American_IMCHIP/release/pc10.txt", sep=" ", header=TRUE, comment.char="~")
fam = read.table("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/AA_icaa_hg19_noIBD.fam", sep=" ", comment.char="~"); names(fam) = c("FAM", "IID", "PAT_ID", "MAT_ID", "sex", "case")
cohort = read.table("/m/jdrfdn_scratch/users/projects/IMCHIP/AA_cohort_list.txt", header=TRUE, comment.char="~")
sources = read.csv("/m/jdrfdn_scratch/users/projects/IMCHIP/T1DGC.2011.03_Golden_Pedigree/T1DGC.2011.03_Golden_Pedigree/Resources/T1DGC.2011.03_Resources.csv")
sources$Race.Code.1[sources$Race.Code.1=="."] <- NA

##Format and merge data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
row.names(pcs) = pcs$IID
row.names(fam) = fam$IID
row.names(cohort) = cohort$ID
row.names(sources) = sources$Analytic.ID

data = data.frame(pcs, fam[row.names(pcs),], cohort[row.names(pcs),])
data$case_fac = factor(data$case, levels=c(1,2), labels=c("control", "case"))
data$COHORT_nomiss = data$COHORT
data$COHORT_nomiss[is.na(data$COHORT) | data$COHORT%in% c("UK", "EUR") ] <- "T1DGC"
data$sex_fac = factor(data$sex, levels=c(1,2), labels=c("females", "males"))

table(data$COHORT_nomiss)
table(data$COHORT_nomiss, data$case)
ids_in_sources_file = data$IID[data$IID%in%sources$Analytic.ID]
table(data[ids_in_sources_file,]$COHORT_nomiss)



##Basic PCA plots - diagnostic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/pca_plots.pdf")
pairs(data[,3:7])
pairs(data[,8:10])


p1 = ggplot(data, aes(x=PC1, y=PC2, col=case_fac)) + geom_point(aes(alpha=0.5)) + geom_point(data=subset(data, case==2), aes(alpha=0.5))+ theme(legend.position="none")
p2 = ggplot(data, aes(x=PC1, y=PC3, col=case_fac)) + geom_point(aes(alpha=0.5)) + geom_point(data=subset(data, case==2), aes(alpha=0.5))+ theme(legend.position="none")
p3 = ggplot(data, aes(x=PC2, y=PC3, col=case_fac)) + geom_point(aes(alpha=0.5)) + geom_point(data=subset(data, case==2), aes(alpha=0.5))+ theme(legend.position="none")
p4 = ggplot(data, aes(x=as.factor(case), y=PC1)) + geom_boxplot()
multiplot(p1, p2, p3, p4, cols=2)


scatterplot3d(data$PC1, data$PC2, data$PC3)
data$case_col = NA
data$case_col[data$case==1] <- "red"
data$case_col[data$case==2] <- "green"
with(data, {scatterplot3d(PC1, PC2, PC3, color=case_col, pch=19)})
color_transparent <- adjustcolor(data$case_col, alpha.f = 0.3)
with(data, {scatterplot3d(PC1, PC2, PC3, color=color_transparent, pch=19)})

#we do see an association between PC1 and case control status, suggesting that cases and controls are not representative
#of exactly the same populations

#try to generate ancestry-matched subset of controls?
#validate finding in subset of cases which can be appropriately ancestry-matched to controls


##PCA plots by source/study/sex/etc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggplot(data, aes(x=PC1, y=PC2, col=COHORT_nomiss=="T1DGC")) + geom_point(aes(alpha=0.5)) #T1DGC vs other cohorts
ggplot(data, aes(x=PC1, y=PC2, col=sex_fac)) + geom_point(aes(alpha=0.5)) + facet_wrap(~COHORT_nomiss, nrow=2) #by cohort 
ggplot(data, aes(x=PC1, y=PC2, col=case_fac)) + geom_point(aes(alpha=0.5)) + facet_wrap(~COHORT_nomiss, nrow=2) #by cohort
ggplot(data, aes(x=PC1, y=PC2, col=case_fac)) + geom_point(aes(alpha=0.5)) + geom_point(data=subset(data, case==2), aes(alpha=0.5)) + facet_wrap(~COHORT_nomiss, nrow=3) #by cohort

td = data[data$COHORT_nomiss=="T1DGC",]
td_source = merge(td, sources, by.x="IID", by.y="Analytic.ID" ); dim(td_source)
td_source$tribal_groups = factor(td_source$Race.Code.1, 
		levels=c("811", "900", "910", "913", "919", "920", "923", "934"),
		labels=c("African-American", "Sub-Saharan","Central_and_West_1", "Ghanaian", "Central_and_West_2", "Southern_and_Eastern", "Eritean", "Somali")
)
table(td_source$tribal_groups)
ggplot(td_source, aes(x=PC1, y=PC2)) +  geom_point(data=td_source[td_source$tribal_groups=="African-American",], color="grey", aes(alpha=0.9, col=tribal_groups)) +
		geom_point(data=td_source[!td_source$tribal_groups=="African-American",], aes(col=tribal_groups)) +
		ggtitle("T1DGC subjects by tribal group")

dev.off()

