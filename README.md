# Master Thesis: Synonymous Mutations in AML

This repository contains all analysis scripts used in my Master's thesis titled:

**"Synonymous mutations in recurrently mutated genes in AML and their impact on Complete Remission and Overall Survival"**

Submitted to **Ludwig Maximilian University of Munich**, Institute for Medical Informatics, Biometry and Epidemiology (IBE).  
Supervisor: Prof. Dr. Ulrich Mansmann  
Reviewer: PD Dr. Markus Pfirrmann  

## 📄 Thesis Abstract
The study investigates the role of synonymous mutations in 68 recurrently mutated genes in Acute Myeloid Leukemia (AML). Using whole-exome sequencing data from the Beat AML cohort, the analysis explores the association between synonymous mutation frequency and patient outcomes such as Complete Remission (CR) and Overall Survival (OS). While no significant association was found with OS, a significant negative association was observed between synonymous mutation burden and CR when stratified by ELN 2017 classification.

## 🧬 Project Overview
This repository includes:
- **R scripts** for statistical analysis, visualization, and data processing
- **AWK scripts** for text processing and variant filtering
- Supplementary analysis pipelines used for variant annotation, filtering, and survival modeling

## 🔬 Methods Overview
1. **Data Source**: Beat AML cohort (175 patients after filtering)
2. **Variant Calling**: VarScan2
3. **Annotation**: ANNOVAR, SynMICdb, CIVIC
4. **Filtering**: GEMINI, bcftools, custom AWK scripts
5. **Statistical Analysis**:
   - Logistic regression for Complete Remission
   - Cox proportional hazards for Overall Survival
   - Stratification by ELN 2017 classification

## 📊 Key Results
- Synonymous mutations were asymmetrically distributed across gene regions.
- Higher synonymous mutation burden was associated with lower Complete Remission rates (OR ≈ 0.60, p ≈ 0.04).
- No significant association was found with Overall Survival.
- Notable mutations were found in `NOTCH1` and `FLT3`.

## 🛠️ Tools & Technologies
- **R**: Statistical analysis, visualization (ggplot2, survival, etc.)
- **AWK**: Text processing for VCF files
- **Bioinformatics Tools**: BWA, GATK, VarScan2, ANNOVAR, GEMINI
- **Databases**: SynMICdb, CIVIC, COSMIC, dbNSFP, 1000 Genomes


## 💬 Discussion & Interpretation:

## 🧬 Key Findings
Asymmetric Distribution: Synonymous mutations were not evenly distributed across gene regions, suggesting functional or structural constraints in mutation tolerance.

Impact on Treatment Response: Higher synonymous mutation burden was associated with significantly lower Complete Remission rates (OR ≈ 0.60, p ≈ 0.04), especially when stratified by ELN 2017 classification.

No Survival Link: No significant association was found with Overall Survival, indicating that synonymous mutations may affect early treatment response rather than long-term survival.

Mutation Hotspots: Clusters were identified in key genes like:

NOTCH1 – signaling and differentiation

U2AF2 – RNA splicing machinery

FLT3 – tyrosine kinase signaling

## ⚙️ Potential Mechanisms
Although “silent” at the protein level, synonymous mutations may influence:

✅ Splicing regulation (exon skipping, splice site disruption)

✅ mRNA stability and translational kinetics

✅ Cotranslational folding and post-translational modification sites

✅ Non-coding RNA binding (e.g., miRNA targeting)

For example:

The NOTCH1 mutation chr9:139395259 is located in an ankyrin repeat domain near known phosphorylation sites, which could alter protein interactions or stability.

## ⚠️ Limitations
Small cohort size limits statistical power.

No independent validation cohort available.

Functional validation (e.g., splicing assays) was not performed.

Rare variant frequency restricted in-depth per-mutation analysis.

## 🚀 Future Directions
🔬 Larger multi-omics cohorts with RNA-seq to assess splicing impacts.

🧪 Functional studies using minigene reporters or CRISPR editing.

📊 Integration with proteomic data to assess protein-level effects.

🏥 Clinical translation into risk stratification tools.

## ✅ Conclusion
This study suggests that synonymous mutations are not biologically silent in AML and may collectively influence therapeutic response. Their inclusion in genomic profiling could refine risk prediction and personalize treatment strategies in AML precision oncology.

