## Generate table for cBioPortal. 
## X-squared test to show that the distribution of variants position are not uniform.  
## Plot the mutation position on their respective genes. 

# Creation of the table for cBioPortal
Sample_ID <- c("09-00705", "14-00081", "14-00787", "09-00705", "14-00175", "14-00739", "13-00273", "15-00610", "15-00680", "14-00832", "16-00001", "16-00056", "14-00329", "14-00464", "14-00725", "16-00003", "16-00035", "14-00096", "15-00839", "16-00001", "14-00632", "14-00781", "15-00813", "09-00705", "10-00172", "10-00507", "10-00692", "13-00165", "13-00195", "13-00253", "13-00273", "13-00425", "13-00466", "13-00513", "13-00659", "14-00034", "14-00175", "14-00253", "14-00259", "14-00423","14-00498", "14-00504", "14-00567", "14-00578", "14-00632", "14-00664", "14-00667", "14-00725", "14-00739", "14-00781", "14-00787", "15-00148", "15-00965", "16-00056", "16-00094", "14-00096", "16-00001", "14-00253", "14-00496", "15-00850", "10-00507", "13-00253", "14-00578", "14-00521", "15-00148", "15-00717", "16-00003", "16-00075", "14-00096", "14-00632", "14-00781", "15-00045", "14-00175", "14-00329", "14-00787", "15-00269", "15-00492", "15-00680", "13-00165", "14-00423", "14-00725", "16-00075", "13-00218", "14-00359", "14-00787", "14-00329", "14-00498", "14-00832", "15-00276", "15-00331", "15-00819", "16-00026", "16-00046", "16-00001", "14-00096", "09-00705", "13-00487", "14-00832", "15-00331", "15-00492", "15-00756", "15-00813", "09-00705", "14-00798", "15-00276", "15-00610", "15-00895", "15-00965", "10-00172", "13-00253", "15-00756", "15-00895", "13-00425", "13-00625", "13-00659", "15-00148", "15-00320", "10-00172", "15-00702", "15-00786", "16-00075", "14-00725", "15-00051", "15-00614", "15-00726", "15-00965", "16-00001", "09-00705", "14-00832", "15-00610", "15-00734", "15-00858", "09-00705", "15-00610", "15-00734", "15-00858", "09-00705", "15-00702", "15-00763", "15-00756", "15-00839", "15-00850", "15-00051", "15-00763", "16-00129")


Cancer_Type <- rep("AML", 145)

Chromo_seq <- c(seq(1:22), "X", "Y")

Chromosome <- c(rep(Chromo_seq[2], 6), rep(Chromo_seq[4], 6), rep(Chromo_seq[5], 5), rep(Chromo_seq[9], 52), rep(Chromo_seq[10], 3), rep(Chromo_seq[11], 10), rep(Chromo_seq[13], 13), rep(Chromo_seq[15], 3), rep(Chromo_seq[18], 4), rep(Chromo_seq[19], 25), rep(Chromo_seq[20], 9), rep(Chromo_seq[21], 3), rep(Chromo_seq[23], 6))

# Start_Position <- c(rep(25457198), rep(25457233, 3), rep(89034550, 3), rep(106190861, 3), rep(1279438, 5), rep(139390675, 3), rep(139391542, 3), rep(139395258, 34), rep(139402714, 3), rep(139440220, 6), rep(8115829, 3), rep(32456725, 6), rep(118374646, 4), rep(28609650, 13), rep(90628536, 3), rep(42643368, 4), rep(33792747, 6), rep(56172452, 4), rep(56180937, 5), rep(56181018, 4), rep(56185319, 6), rep(31023499, 5), rep(31024033, 4), rep(36164485, 3), rep(129148811, 3))

Start_Position <- End_Position

End_Position <- c(rep(25457199, 3), rep(25457234, 3), rep(89034551, 3), rep(106190862, 3), rep(1279439, 5), rep(139390676, 3), rep(139391543, 3), rep(139395259, 34),rep(139397728, 3),  rep(139402715, 3), rep(139440221, 6), rep(8115830, 3), rep(32456726, 6), rep(118374647, 4), rep(28609651, 13), rep(90628537, 3), rep(42643369, 4), rep(33792748, 6), rep(56172453, 4), rep(56180938, 5), rep(56181019, 4), rep(56185320, 6), rep(31023500, 5), rep(31024034, 4), rep(36164486, 3), rep(15838411, 3), rep(129148812, 3))



Reference_Allele <- c(rep("T", 3), rep("T", 3), rep("C", 3), rep("T", 3), rep("G", 5), rep("A", 3), rep("C", 3), rep("G", 34), rep("C", 3), rep("G", 3), rep("C", 6), rep("C", 3), rep("G", 6), rep("C", 4), rep("C", 13), rep("G", 3), rep("C", 4), rep("G", 6), rep("C", 4), rep("C", 5), rep("C", 4), rep("T", 6), rep("C", 5), rep("G", 4), rep("G", 3), rep("C", 3), rep("C", 3))


Variant_Allele <- c(rep("C", 3), rep("G", 3), rep("T", 3), rep("C", 3), rep("A", 5), rep("C", 3), rep("T", 3), rep("A", 34), rep("T", 3), rep("A", 3), rep("A", 6), rep("T", 3), rep("T", 6), rep("T", 4), rep("G", 13), rep("A", 3), rep("T", 4), rep("A", 6), rep("T", 4), rep("T", 5), rep("T", 4), rep("A", 6), rep("T", 5), rep("A", 4), rep("C", 3), rep("T", 3), rep("T", 3))



Hugo_Symbol <- c(rep("DNMT3A", 6), rep("ABCG2", 3), rep("TET2", 3), rep("TERT", 5), rep("NOTCH1", 52), rep("GATA3", 3), rep("WT1", 6), rep("KMT2A", 4), rep("FLT3", 13), rep("IDH2", 3), rep("SETBP1", 4), rep("CEBPA", 6), rep("U2AF2", 4), rep("U2AF2", 5), rep("U2AF2", 4), rep("U2AF2", 6), rep("ASXL1", 5), rep("ASXL1", 4), rep("RUNX1", 3), rep("ZRSR2", 3), rep("BCORL1", 3))



Protein_Change <- c(rep("P896P", 3), rep("R885R", 3), rep("E366E", 3), rep("H1401H", 3), rep("A699A", 5), rep("P2505P", 3), rep("P2216P", 3), rep("G1893G", 34), rep("Q1691Q", 3), rep("S1098S", 3), rep("A6A", 6), rep("N393N", 3), rep("R56R", 6), rep("D2680D", 4), rep("T526T", 13), rep("T350T", 3), rep("A1499A", 4), rep("H191H", 6), rep("D128D", 4), rep("D391D", 5), rep("P418P", 4), rep("S438S", 6), rep("H995H", 5), rep("L1173L", 4), rep("P463P", 3), rep("P303P", 3), rep("V688V", 3))


Mutation_Type <- rep("Synonymous_Mutation", 145)


genomic_location <- cbind(Sample_ID, Cancer_Type,Hugo_Symbol, Protein_Change, Mutation_Type, Chromosome, Start_Position, End_Position, Reference_Allele, Variant_Allele)

write.table(genomic_location, 'genom_locat.txt', col.names  = T, row.names = F, quote = F, sep = "\t")

C_2_DNMT3A_1 <- 896/912
C_2_DNMT3A_2 <- 885/912 
C_4_ABCG2 <- 366/615 
C_4_TET2 <- 1401/2023 
C_5_TERT <- 699/1132 
C_9_NOTCH1_1 <- 2505/2555 
C_9_NOTCH1_2 <- 2216/2555 
C_9_NOTCH1_3 <- 1893/2555 
C_9_NOTCH1_4 <- 1691/2555
C_9_NOTCH1_5 <- 1098/2555 
C_9_NOTCH1_6 <- 6/2555 
C_10_GATA3 <- 393/444 
C_11_WT1 <- 56/571 
C_11_KMT2A <- 2680/3972 
C_13_FLT3 <- 526/996 
C_15_IDH2 <- 350/452 
C_18_SETBP1 <- 1499/1596 
C_19_CEBPA <- 191/358 
C_19_U2AF2_1 <- 128/475 
C_19_U2AF2_2 <- 391/475 
C_19_U2AF2_3 <- 418/475 
C_19_U2AF2_4 <- 438/475 
C_20_ASXL1_1 <- 995/1541 
C_20_ASXL1_2 <- 1173/1541 
C_21_RUNX1 <- 463/480
C_X_ZRSR2 <- 303/482
C_X_BCORL1 <- 688/1785 

reScale_Values <- c(rep(C_2_DNMT3A_1, 3), rep(C_2_DNMT3A_2, 3), rep(C_4_ABCG2, 3), rep(C_4_TET2, 3), rep(C_5_TERT, 5), rep(C_9_NOTCH1_1, 3), rep(C_9_NOTCH1_2, 3), rep(C_9_NOTCH1_3, 34), rep(C_9_NOTCH1_4, 3), rep(C_9_NOTCH1_5, 3), rep(C_9_NOTCH1_6, 6), rep(C_10_GATA3, 3), rep(C_11_WT1, 6), rep(C_11_KMT2A, 4), rep(C_13_FLT3, 13), rep(C_15_IDH2, 3), rep(C_18_SETBP1, 4), rep(C_19_CEBPA, 6), rep(C_19_U2AF2_1, 4), rep(C_19_U2AF2_2, 5), rep(C_19_U2AF2_3, 4), rep(C_19_U2AF2_4, 6), rep(C_20_ASXL1_1, 5), rep(C_20_ASXL1_2, 4), rep(C_21_RUNX1, 3), rep(C_X_ZRSR2, 3), rep(C_X_BCORL1, 3))


densAna <- data.frame(Sample_ID, Hugo_Symbol, Protein_Change, reScale_Values, Chromosome, End_Position, Reference_Allele, Variant_Allele)

# saving the file.
save(file = "Density.RData", densAna)
setwd("C:/AAA LMU/4th/Thesis/Analysis")
load(file = "Density.RData")

hist(unique(densAna$reScale_Values))
hist(unique(densAna$reScale_Values), 
     main = "Relative distribution of variants position in protein sequence", 
     xlab = "Variants position in protein sequence", 
     ylab = "Number of variants", 
     border = "blue", 
     col = "red", breaks = 20)

library(dplyr)
library(ggplot2)

tmp2 <- densAna %>% group_by(Protein_Change, reScale_Values) %>% count() %>% ungroup()
## Here rescale values are unique

qplot(reScale_Values, data = tmp2, geom = "histogram", binwidth = 0.02, fill=I("blue"), col=I("red"), alpha=I(.7), xlab = "Relative position of the mutations", ylab = "Number of mutations", main = "Distribution of mutations position in protein sequence (Unique)")
# lines(density(tmp2$reScale_Values))


# 
# qplot(tmp2$reScale_Values, stat = "ecdf", geom = "step", xlab = "Variants position", ylab = "Distribution function of variant position", main = "Empirical cumulative distribution function of variants position")
# ## (** stat is deprecated) 

plot(ecdf(tmp2$reScale_Values), xlab="Variants position in protein sequence", ylab = "Distribution function f(x)", do.points = FALSE, main="Empirical cumulative distribution function of variants position")
## The distribution is clearly left skewed.


## X-squared test:
cont_table <- data.frame(c(13.5, 5), c(13.5, 22))
chisq.test(cont_table)


library(ggplot2)

qplot(reScale_Values, data = densAna, geom = "histogram", binwidth = 0.02, fill=I("blue"), col=I("red"), alpha=I(.7), xlab = "Relative position of the mutations", ylab = "Number of mutations", main = "Distribution of mutations position in protein sequence (Cumulative)")
#lines(density(densAna$reScale_Values))

# 
# plot(unique(densAna$reScale_Values))
# adb <- unique(densAna$reScale_Values)
# plot(adb, pch = 19, col = "black")
