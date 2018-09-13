options(stringsAsFactors=FALSE)
setwd("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno/hla_assoc")


args=commandArgs(trailingOnly=TRUE)
if (length(args)==2) {
	genofile = args[1]
	pdffile = args[2]
	#genofile = "derived_data/genotypes_recode.raw.2d"
	#pdffile = "tables/hwe_plots.2d.pdf"
} else {
	stop("Incorrect number of arguments\n")
}

#genofile = "derived_data/genotypes_recode.raw.2d"
#pdffile = "tables/hwe_plots.2d.pdf"
geno = read.table(genofile, comment.char="~", header=TRUE)
genes = c("HLA_A", "HLA_B", "HLA_C", "HLA_DRB1", "HLA_DQA1", "HLA_DQB1", "HLA_DPA1", "HLA_DPB1")


#Hardy-Weinberg test
getObsExpTable = function(x) {
	n= dim(x)[1]
	p = NULL
	for (i in 1:length(alleles)) {
		p[i] = sum(x[,alleles[i]])/(2*n)
	}
	
	expected = NULL
	allele1 = NULL
	allele2 = NULL
	p1_ = NULL
	p2_ = NULL
	count=1
	observed = NULL
	for (i in 1:length(alleles)) {
		p1 = p[i]
		for (j in 1:length(alleles)) {
			p2=p[j]
			if (i==j) {
				e = (p1^2)*n
				allele1[count] = alleles[i]
				allele2[count] = alleles[j]
				expected[count] = e
				observed[count] = sum(x[,allele1[count]]==2)
				p1_[count] = p1
				p2_[count] = p2
				count = count+1
			} else if (i<j) {
				e = 2*p1*p2*n
				allele1[count] = alleles[i]
				allele2[count] = alleles[j]
				expected[count] = e
				observed[count] = sum(x[,allele1[count]]==1 & x[,allele2[count]]==1)
				p1_[count] = p1
				p2_[count] = p2
				count = count+1
			}
		}
	}
	data.frame(allele1, p1_, allele2, p2_, expected, observed)
}

testHW = function(tbl, c) {
	chi = NULL
	#c = 0.5
	m = dim(tbl)[1]
	for (i in 1:m) {
		exp = tbl[i,"expected"]
		obs = tbl[i,"observed"]
		chi[i] = (abs(obs-exp)-c)^2/exp
	}
	pchisq(sum(chi), df=m*(m-1)/2)
}



pdf(pdffile)
for (i in 1:length(genes)) {
	gene = genes[i]
	alleles=names(geno)[grep(gene, names(geno))]
	X = geno[,alleles]
	y = geno$PHENOTYPE-1 #case=2 control=1
	
	#Quality filtering
	#summary(rowSums(X))
	X_filt = X[rowSums(X)==2,]
	y_filt = y[rowSums(X)==2]
	
	#HWE tests #alternatives to chisq test? http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1415-47572009000300028
	case_tbl <- getObsExpTable(X_filt[y_filt==1,])
	control_tbl <- getObsExpTable(X_filt[y_filt==0,])
	common = control_tbl[control_tbl[,"p1_"]>0.01&control_tbl[,"p2_"]>0.01,]
	
	par(mfrow=c(2,2))
	plot(control_tbl[,"expected"], control_tbl[,"observed"], main=paste(gene,"\nControls, p=", 1-testHW(control_tbl, c=0.5), sep=""))
	plot(case_tbl[,"expected"], case_tbl[,"observed"], main=paste("Cases, p=", 1-testHW(case_tbl, c=0.5), sep=""))
	plot(control_tbl[,"observed"], case_tbl[,"observed"], main="Cases vs Controls")
	plot(common[,"expected"], common[,"observed"], main=paste("Common alleles only, p=", 1-testHW(common, c=0.5), sep=""))
}
dev.off()