options(stringsAsFactors=FALSE)
setwd("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_assoc")


runConditionalAnalysis = function(infile, threshold) {

	dose = read.table(infile)
	fam = read.table("../hla_imputation/imputation_results/AA_IMPUTED.fam", comment.char="~")
	names(fam) <- c("FID","IID", "FAT", "MOT", "SEX", "CASE") #case=2 control=1
	uniq_ids = paste(fam$FID, fam$IID, sep="_")
		
	getGLMStats = function(glm_obj, ind_var) {
		coefs = coef(summary(glm_obj))
		est = coefs[ind_var,"Estimate"] 
		or = exp(est)
		p = coefs[ind_var,4]
		c(est, or, p)
	}
	
	
	runAssoc = function(x, y) {	
		out_glm = glm(y~x, family="binomial")
		case_freq = mean(x[y==1])
		control_freq = mean(x[y==0])
		freq = mean(x)
		c(case_freq, control_freq, freq, getGLMStats(out_glm, "x"))
	}
	
	
	runCondAssoc = function(x, y, condition_mat) {	
		X_new = cbind(x, condition_mat)
		out_glm = glm(y~X_new, family="binomial")
		case_freq = mean(x[y==1])
		control_freq = mean(x[y==0])
		freq = mean(x)
		c(case_freq, control_freq, freq, getGLMStats(out_glm, "X_newx"))
	}
			
	getInfo = function(assoc_mat) {
		min_p = min(assoc_mat[,"P"], na.rm=TRUE)
		topVar = rownames(assoc_mat)[assoc_mat[,"P"]==min_p]
		list(topVar, assoc_mat[topVar,])
	}
	
	#Define input matrices
	X <- t(as.matrix(dose[,4:dim(dose)[2]])); colnames(X) <- dose[,1]; rownames(X) <- uniq_ids
	y <- fam$CASE-1
	
	#Initialize loop	
	min_p=0
	condition_list_all=c()
	condition_list_pruned=c()
	output = list()
	
	round=1
	while (min_p<threshold) {
		cat(round,"\n")
		
		#run model
		if (length((condition_list_all))==0) {
			out = t(apply(X=X, FUN=runAssoc, MARGIN=2, y))   #<<<<---- INSERT COVARIATES INTO X matrix
		} else {
			out = t(apply(X=X[,!colnames(X)%in%condition_list_all], FUN=runCondAssoc, MARGIN=2, y, X[,condition_list_pruned]))
		}
		colnames(out) <- c("case_AF", "control_AF", "total_AF", "Est","OR", "P")
		#export association results for each round
		output[[round]] = list(round=round, condition_snps=condition_list_all, strongest_variant=getInfo(out)[[1]], assoc_stats=getInfo(out)[[2]], out)
		
		#update		
		if (class(output[[round]]$assoc_stats)=="numeric"){
			min_p = output[[round]]$assoc_stats["P"]	
		} else if (class(output[[round]]$assoc_stats)=="matrix") {
			min_p = output[[round]]$assoc_stats[,"P"][1]
		}
		condition_list_all = c(condition_list_all, output[[round]]$strongest_variant)
		condition_list_pruned = c(condition_list_pruned, output[[round]]$strongest_variant[1])
		round = round+1		
	}
	output
}

getSummaryTable = function(output, outfile) {
	print_out = data.frame(unlist(lapply(output, '[[', 3)), matrix(unlist(lapply(output, '[[', 4) ),ncol=6,byrow=TRUE))
	names(print_out) <- c("VAR", "case_AF", "control_AF", "total_AF", "Est","OR", "P")
	write.table(print_out, file=outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
	print(print_out)
}

out.2d =runConditionalAnalysis(infile = "processed_impdata/dosages.2d", threshold=5e-8)
getSummaryTable(out.2d, outfile="assoc/conditional_analysis_2d.out")

out.4d =runConditionalAnalysis(infile = "processed_impdata/dosages.4d", threshold=5e-8)
getSummaryTable(out.4d, outfile="assoc/conditional_analysis_4d.out")
sumtable = getSummaryTable(out.aa, outfile="assoc/conditional_analysis_aa.out")
sumtable$gene = sapply(strsplit(sumtable$VAR, split="_"),"[[", 2)
#get amino acid changes meeting bonferonni significance
sumtable[sumtable$P<3e-4,] 
##Compare these sites with EUR
sites=m_short[m_short$SNP%in%sumtable$VAR,]
sites[sites$P_EUR>-3,]

#out.aa =runConditionalAnalysis(infile = "processed_impdata/dosages.aa", threshold=0.05)
#save(out.aa, file="assoc/conditional_analysis_object.aa.RData")
load(file="conditional_analysis_object.aa.RData")
getSummaryTable(out.aa, outfile="assoc/conditional_analysis_aa.out")

out.all =runConditionalAnalysis(infile = "processed_impdata/dosages.all", threshold=5e-5)
getSummaryTable(out.all, outfile="assoc/conditional_analysis_all.out")



### Plot Amino Acid associations
library(ggplot2, lib.loc="/home/ccr5ju/R/x86_64-redhat-linux-gnu-library/3.3")




getData = function(d) {
	varsplit = strsplit(row.names(d), split="_")
	d$gene = sapply(varsplit,"[[", 2)
	d$Gene = factor(d$gene, levels=c("A","C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"), labels=c("HLA-A","HLA-C", "HLA-B", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1"))
	d$aa_pos = sapply(varsplit,"[[", 3)
	d$aa_pos_numeric=d$aa_pos
	d$aa_pos_numeric[grep("x", d$aa_pos)] = sapply(strsplit(d$aa_pos[grep("x", d$aa_pos)],split="x"),"[[",1)
	d$aa_pos_numeric = as.numeric(d$aa_pos_numeric)
	#d$bp_pos = as.numeric(sapply(varsplit,"[[", 4))
	d
}


##### PLOTS FOR POSTER
makePlots = function(df, title, xlab) {
	axisstyle.y <- element_text(face = "bold", size = 15)
	axisstyle.x <- element_text(size = 12)
	axislabelstyle.x =element_text(size=16)
	ggplot(df, aes(x=aa_pos_numeric, y=-log10(P))) + geom_point() + facet_grid(.~Gene) + 
			ggtitle(title) + ylab("-log10(p)") + xlab(xlab)+ #xlab("Amino acid position")
			theme(axis.text.y = axisstyle.y, axis.text.x = axisstyle.x, axis.title.x=axislabelstyle.x) +
			theme(aspect.ratio=1/2) + 
			theme(strip.text.x = element_text(size = 12, face="bold"))
}

### Round 1 plot
pdf("assoc/amino_acid_conditional_analysis_round1.pdf")
	d = getData(as.data.frame(out.aa[[1]][[5]]))
	makePlots(d[d$gene%in%c("DRB1", "DQA1", "DQB1"),], title="", xlab="")
dev.off()

### Round 2 plot
pdf("assoc/amino_acid_conditional_analysis_round2.pdf")
	d = getData(as.data.frame(out.aa[[2]][[5]]))
	makePlots(d[d$gene%in%c("DRB1", "DQA1", "DQB1"),], title="",xlab="")
dev.off()

### Round 3 plot
pdf("assoc/amino_acid_conditional_analysis_round3.pdf")
	d = getData(as.data.frame(out.aa[[3]][[5]]))
	makePlots(d[d$gene%in%c("DRB1", "DQA1", "DQB1"),], title="", xlab="Amino acid position")
dev.off()

#scp assoc/amino_acid_conditional_analysis_round*.pdf "$local/AA_Immunochip"



##### PLOTS FOR EXPLORATORY ANALYSIS
makePlots2 = function(df, title) {
	axisstyle.y <- element_text(face = "bold", size = 8)
	axisstyle.x <- element_text(size = 8)
	ggplot(df, aes(x=aa_pos_numeric, y=-log10(P))) + geom_point() + facet_grid(.~Gene) + 
			ggtitle(title) + ylab("-log10(p)") + xlab("")+ #xlab("Amino acid position")
			theme(axis.text.y = axisstyle.y, axis.text.x = axisstyle.x) +
			theme(aspect.ratio=1/2) + 
			theme(strip.text.x = element_text(size = 12, face="bold"))
}

#All rounds, all genes
pdf("assoc/amino_acid_conditional_analysis_all.pdf")
for (i in 1:length(out.aa)) {
#for (i in 1:5) {
	d = getData(as.data.frame(out.aa[[i]][[5]]))
	p = makePlots2(d, title="")
	print(p)
}
dev.off()
#scp assoc/amino_acid_conditional_analysis_all.pdf "$local/AA_Immunochip"


#associations with remaining variants
dose = read.table("processed_impdata/dosages.all")
aa_condition_list = sumtable$VAR[!sumtable$VAR%in%c("AA_DQA1_52_32717206_H","AA_DQA1_54_32717212")]

fam = read.table("../hla_imputation/imputation_results/AA_IMPUTED.fam", comment.char="~")
names(fam) <- c("FID","IID", "FAT", "MOT", "SEX", "CASE") #case=2 control=1
uniq_ids = paste(fam$FID, fam$IID, sep="_")

X <- t(as.matrix(dose[,4:dim(dose)[2]])); colnames(X) <- dose[,1]; rownames(X) <- uniq_ids
y <- fam$CASE-1

runCondAssoc = function(x, y, condition_mat) {	
	if (max(cor(x, condition_mat))<0.95) {
		X_new = cbind(x, condition_mat)
		out_glm = glm(y~X_new, family="binomial")
		case_freq = mean(x[y==1])
		control_freq = mean(x[y==0])
		freq = mean(x)
		out=c(case_freq, control_freq, freq, getGLMStats(out_glm, "X_newx"))
	} else {
		case_freq = mean(x[y==1])
		control_freq = mean(x[y==0])
		freq = mean(x)
		out=c(case_freq, control_freq, freq, NA, NA, NA)
	}
	
}
getGLMStats = function(glm_obj, ind_var) {
	coefs = coef(summary(glm_obj))
	est = coefs[ind_var,"Estimate"] 
	or = exp(est)
	p = coefs[ind_var,4]
	c(est, or, p)
}
out = as.data.frame(t(apply(X=X[,!colnames(X)%in%aa_condition_list], FUN=runCondAssoc, MARGIN=2, y, X[,aa_condition_list])))
colnames(out) <- c("case_AF", "control_AF", "total_AF", "Est","OR", "P")
out$VAR = row.names(out)
write.table(out, "assoc/conditioned_on_all_aa.txt")
#sed 's/"//g' conditioned_on_all_aa.txt | awk 'NR==1 {print "SNP",$0} NR>1 {print}' > conditioned_on_all_aa_errorsRemoved.txt
#scp conditioned_on_all_aa_errorsRemoved.txt "$local/AA_Immunochip"

out_sorted = out[order(out$P),]


((0.4141)/(1-0.4141))/((0.5507)/(1-0.5507))
