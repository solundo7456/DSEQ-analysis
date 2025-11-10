# DSEQ-analysis
This project demonstrates a complete RNA-seq differential expression analysis pipeline using publicly available SARS-CoV-2 transcriptomic data (GEO: GSE152075). The goal is to identify genes that are differentially expressed between SARS-CoV-2 positive and negative samples.

Key steps in the analysis include:

Data Loading and Inspection

Raw counts loaded from GEO and inspected for quality.

Metadata extracted and cleaned for downstream analysis.

Quality Control & Filtering

Library sizes assessed and visualized.

Lowly expressed genes filtered to improve statistical power.

Data Preparation

Metadata and count matrices aligned.

Factors and numeric variables properly formatted for analysis.

Differential Expression Analysis

DESeq2 used to model gene expression differences based on infection status.

Significantly differentially expressed genes identified (FDR < 0.05).

Visualization

MA-plot for overview of gene expression changes.

Volcano plot highlighting significant DE genes.

Heatmap of top 50 DE genes across samples.

Results Export

Differential expression results exported for downstream exploration and reproducibility.

This project is ideal for bioinformaticians, data scientists, and researchers interested in transcriptomic analysis, viral infection studies, and reproducible RNA-seq workflows.

Tools & Libraries: R, DESeq2, limma, GEOquery, EnhancedVolcano, pheatmap, statistical modeling, data wrangling.
