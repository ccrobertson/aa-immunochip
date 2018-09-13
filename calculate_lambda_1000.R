library(qqman)
setwd("/m/jdrfdn_scratch/users/ccr5ju/AA_Immuno")
datafile = "/m/jdrfdn_scratch/users/so4g/IMCHIP/20160316_African_American_IMCHIP/imputation/association/AA_immuno_imputed_merged.assoc.logistic"
d = read.table(datafile, header=TRUE)

getPlots = function(d, title) {
	CHISQ <- qchisq(1-d$P,1)
	lambda = median(CHISQ)/qchisq(0.5,1)
	ncase = 1036
	ncontrol = 2313
	lambda_1k = 1 + (lambda-1)*((1/ncase + 1/ncontrol)/(1/1000 + 1/1000))
	qq(d$P, main = paste(title,"\n Lambda_1000=", round(lambda_1k, digits=2)))
	lambda_1k
}

pdf("qq_plots.pdf")
getPlots(d, title="")

d_trimmed_1 = d[!((d$CHR==6 & d$BP>20000000 & d$BP<40000000)|(d$CHR==11 & d$BP>1000000 & d$BP<3000000)|(d$CHR==1 & d$BP>112000000 & d$BP<115000000)),] #remove MHC, INS, and PTPN22
getPlots(d_trimmed_1, title="Excluding MHC, INS, and PTPN22 regions")

d_trimmed_2 = d[!((d$CHR==6 & d$BP>20000000 & d$BP<40000000)|(d$CHR==11 & d$BP>1000000 & d$BP<3000000)),] #remove MHC, and INS
getPlots(d_trimmed_2, title="Excluding MHC and INS regions")

dev.off()

#scp qq_plots.pdf "$local/AA_Immunochip/tables_and_figures" 
