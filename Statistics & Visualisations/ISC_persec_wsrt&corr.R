#install.packages("R.matlab")

p_vals_wsrt <- c()
corr_vals <- c()

for (movie in 1:20)
{
  for (comp in 1:3)
  {
    path_drug1 <- paste("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/ISC_CALC/DRUG1/VIDEO",movie,"_ISC_persecond.mat", sep="")
    data_drug1 <- readMat(path_drug1)
    path_drug2 <- paste("/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/ISC_CALC/DRUG2/VIDEO",movie,"_ISC_persecond.mat", sep="")
    data_drug2 <- readMat(path_drug2)
    
    # perform Wilcoxon  signed-rank test
    # paired has to be TRUE; if FALSE it's a Mann-Whitney Test
    res <- wilcox.test(data_drug1[["ISC.persecond.drug1"]][comp,],data_drug2[["ISC.persecond.drug2"]][comp,], paired = TRUE, alternative = "two.sided")
    p_vals_wsrt <- c(p_vals_wsrt, res$p.value)
    
    # perform Pearson Correlation between both time series of drug1 and drug2
    res <- cor(data_drug1[["ISC.persecond.drug1"]][comp,], data_drug2[["ISC.persecond.drug2"]][comp,])
    corr_vals <- c(corr_vals, res)
  }
}

write(p_vals_wsrt, "/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/ISC_CALC/stats/pvals_WSRT.txt",ncolumns = 3, sep = "\t")
write(corr_vals, "/Volumes/LaCie Armaz/DPrataLab/ISC/FINAL_PREPROCESSING/ISC_CALC/stats/corr_vals.txt",ncolumns = 3, sep = "\t")
