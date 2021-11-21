## Data preparation 

## Loading necessary libraries
library(tidyr)
library(dplyr)

# Loading the 'variant calling file (vcf)'. 
raw_data_tmp <- data.table::fread("variants.vcf", header = T, sep = '\t')

# Loading the 46480 synonymous variants, accuired by 'Gemini query' tools.
gs_list <- read.table("filt258.tabular", header = T)

# The 68 recurrently mutated genes in AML. 
genes <- c("ADA", "ABCB1", "ABCG2", "ASXL1", "BCOR", "BCORL1",
           "BRAF", "CBL", "CDA", "CDKN2A", "CEBPA", "CSF3R",
           "CSFR1", "DAXX", "DCK", "DCLK1", "DIS3", "DNMT3A",
           "ETV6", "EZH2", "FAM5C", "FBXW7", "FLT3", "GATA1",
           "GATA2", "GATA3", "HNRNPK", "HRAS", "IDH1", "IDH2",
           "IL7R", "JAK1", "JAK2", "JAK3", "KDM6A", "KMT2A",
           "KIT", "KRAS", "MIR-142", "MPL", "MYD88", "NOTCH1",
           "NPM1", "NRAS", "NT5C2", "PHF6", "PTEN", "PTPN11",
           "PTPRT", "RAD21", "RUNX1", "SETBP1", "SF1", "SF3A1", "SF3B1",
           "SMC1A", "SMC3", "SRSF2", "STAG2", "TERC", "TERT", "TET2", 
           "TP53", "U2AF1", "U2AF2", "WAC", "WT1", "ZRSR2")

# Filter the synonymous variants found only in the 68 genes. 
newGs2 <- gs_list %>% filter(gene %in% genes)
newGs2 <- newGs2 %>% mutate(chrpos = paste(chrom, end, ref, alt, sep = '_'))


gs_list <- as.data.frame(gs_list) %>% mutate(chrpos = paste(chrom, end, ref, alt, sep = '_'))

# Merge the synonymous variants with their variant allele frequency table
raw_data_tmp2 <- raw_data_tmp %>% filter(chrpos %in% gs_list$chrpos)

# Define the allele frequency in three categories: minor, medium, and major
raw_data_tmp2 <- raw_data_tmp2 %>% mutate_at(2:ncol(.), function(e){
  ifelse(e < 2, 0, ifelse(e < 75, 1, 2))
})


# Separate the synonymous mutations belong to the 68 recurrently mutated genes in AML.
# All together 190  have found.
tmp <- raw_data_tmp2 %>% filter(chrpos %in% newGs2$chrpos) %>% as.data.frame()



tmp <- tmp %>% pivot_longer(cols = starts_with("Sample"))
# Set the variable 'in_68genes' 
tmp1 <- tmp %>% group_by(name) %>% summarize(syncount = sum(value, na.rm = T)) %>% mutate(in_68genes = ifelse(syncount > 0, 1, 0), 
                                                                                          sample_id = gsub('_AML.realign.bam', '', gsub('Sample_', '', name)))


## Only to save the data.
write.table(tmp1, 'synonymous_variant_count.txt', row.names = F, col.names = T, quote = F, sep = '\t')

## Load the table of clinical characteristics of the patients; 
c_samples <- read.csv("Clinical_Summary.csv", sep = ";") 

# Filter the patients with the selected criteria
abc <- c_samples %>% filter(specimenGroups == "Initial Acute Leukemia Diagnosis" &
                              isRelapse == "FALSCH" &
                              specificDxAtInclusion != "Acute promyelocytic leukaemia with t(15;17)(q22;q12); PML-RARA" &
                              specificDxAtInclusion != "Myeloid sarcoma" &
                              inferred_ethnicity == "White")

# Exclude double patients ID
cl_sam5 <- abc %>% group_by(patientId) %>% filter(! n()>1)
names(cl_sam5)[1] <- "Sample_id"

tmp1 <- tmp1 %>% filter(Sample_id %in% cl_sam5$Sample_id)
cl_sam5 <- cl_sam5 %>% filter(Sample_id %in% tmp1$Sample_id)

## Merging variant table with the clinical characteristics of the patients.
data <- merge(tmp1, cl_sam5, by = "Sample_id")

## Creating the 'Vital' variable for survival analysis
data <- data %>% mutate(
  dead = ifelse(
    !vitalStatus %in% c("Alive", "Dead"), 
    NA,
    ifelse(vitalStatus == "Alive", 0, 1)
  )
)

# The variable 'syncount' created before is wrong. Thus, delete it and create correctly. 
data <- data %>% select(-syncount)

## New 'syncount' variable-
data_tmp <- data %>% select(grep('^chr', names(data),
                                 value = T)) %>% summarize(syncount = rowSums(. > 0, na.rm = T),
                                                                                syndosage = rowSums(.,
                                                                                                    na.rm = T))%>% mutate(syncount_binary = ifelse(syncount > 0, 1,
                                                                                                                                                   0))%>% select(syncount,
                                                                                                                                                   syncount_binary, syndosage)

data <- cbind(data, data_tmp)

## Crete ELN:2017 classification model stratified.
data <- data  %>% mutate(ELN_2017 = factor(ELN2017, c("Favorable", "Intermediate", "Adverse")))

save(file = "FinalData.RData", data) 

