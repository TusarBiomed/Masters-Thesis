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

- ## 🔬 Methods Overview
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


Key Findings and Interpretation:
Our exploratory analysis revealed that synonymous mutations in recurrently mutated AML genes are not uniformly distributed across protein sequences, suggesting potential functional constraints or region-specific mutational pressures. While individual synonymous mutations did not significantly impact overall survival, the collective burden of such mutations was significantly associated with reduced likelihood of complete remission, particularly when stratified by ELN 2017 risk classification. This indicates that synonymous mutations may exert a cumulative biological effect, possibly through mechanisms such as altered mRNA splicing, stability, or translational efficiency, which collectively influence therapeutic response.

Biological and Clinical Implications:
The presence of synonymous mutation clusters in genes such as NOTCH1 and U2AF2 — both implicated in splicing regulation and leukemogenesis — underscores the potential role of “silent” mutations in RNA processing and protein function. For instance, the recurrent synonymous mutation in NOTCH1 (chr9:139395259) lies within an ankyrin repeat domain near known post-translational modification sites, suggesting possible effects on protein interaction or stability. Similarly, mutations in U2AF2 within RNA recognition motifs may influence spliceosome assembly and fidelity. These findings align with growing evidence that synonymous mutations can modulate gene expression and protein conformation, thereby contributing to disease phenotype and therapy resistance.

Study Limitations and Future Directions:
The primary limitations of this study include its modest sample size, lack of validation cohort, and the rarity of individual synonymous mutations, which constrained robust statistical modeling. Future studies should prioritize larger, multi-center cohorts with matched RNA-seq data to directly assess the impact of synonymous mutations on splicing patterns, transcript stability, and protein expression. Functional assays, such as minigene splicing reporters or CRISPR-based mutagenesis, would be valuable to validate the mechanistic effects of prioritized mutations. Additionally, integrating synonymous mutation profiles with transcriptomic and proteomic datasets may uncover novel biomarkers for risk stratification or therapeutic targeting in AML.

Conclusion:
In conclusion, while synonymous mutations have historically been overlooked in AML pathogenesis, our findings suggest they may collectively influence disease biology and treatment response. This study provides a foundational framework for further investigation into the functional relevance of synonymous mutations in AML and highlights the importance of incorporating non-coding and silent variants into comprehensive genomic profiling for precision oncology.


