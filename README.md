#IBM Comprehensive Transcriptomics Analysis

This repository contains data and code for the data analysis performed in the study: "Comprehensive transcriptomic analysis shows disturbed calcium homeostasis and deregulation of T lymphocyte apoptosis in inclusion body myositis" (Johari et al., 2022, J Neurol. 2022 Aug;269(8):4161-4173. doi: 10.1007/s00415-022-11029-7. Epub 2022 Mar 2. PMID: 35237874; PMCID: PMC9293871).

##Data Analysis
###DESeq2 Analysis
    Main Dataset: 39_PE_samples_fragmentcounts_noMulti_noOverlap_updated.txt
    lncRNA Dataset: lncRNA_39_PE_samples_fragmentcounts_noMulti_noOverlap_updated.txt
    Script: Deseq_analysis_script_johari_etal2021.r

###Venn Diagram (Figure 1)
    Script: Fig1.R

###Gene Expression Plots (Figures 2 & 5)
    Script: Fig2_5.R

###UpSet Plot and ClusterProfiler Analysis (Figure 4)
    Script: Fig4.R

**To replicate the analysis:**

Clone the repository:
```bash
    git clone https://github.com/mriduljohari/IBM-transcriptomics
```    
Navigate to the cloned repository directory and source the R scripts in R within the git repository directory.
```R
    source ("Deseq_analysis_script_johari_etal2021.r");
    source ("Fig1.R");
    source ("Fig2_5.R");
    source ("Fig4.R");
```
For a detailed description of the methods, please refer to the original publication.

This work is part of a research study conducted by Mridul Johari and collaborators. For further information or queries, please contact mridul.johari@helsinki.fi
