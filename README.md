Description
This repository provides an end-to-end workflow for bulk RNA-seq analysis in mouse models, focusing on differential gene expression and functional interpretation. The pipeline includes:
Data Preprocessing
Convert raw .tab files to .csv using Python.
Summarize gene counts and clean metadata in R.
Exploratory Analysis
Perform PCA and UMAP to assess sample distribution.
Identify and remove outlier samples to improve clustering.
Differential Expression Analysis
Conduct DESeq2-based DEG analysis between experimental groups.
Export ordered results for downstream interpretation.
Visualization
Generate customized volcano plots highlighting significant DEGs and GO-term–associated genes.
Create heatmaps of GO-related genes across experimental groups.
Functional Enrichment
Annotate DEGs with GO and KEGG terms.
Highlight biological processes of interest such as chemotaxis, endothelial cell proliferation, and VEGF production.
The workflow is fully annotated, modular, and export-ready (PDF/CSV outputs). It is suitable for reproducible RNA-seq analysis and can be easily adapted to other datasets.
