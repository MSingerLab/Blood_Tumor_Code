Blood_Tumor_Code is the GitHub repository that holds the code for the paper published at the Journal of Experimental Medicine:
"Single-cell analyses identify circulating anti-tumor CD8 T cells and markers for their enrichment."

Find the paper here!
https://doi.org/10.1084/jem.20200920

This directory has four subdirectories:

- DataProcessing
Used for converting 10X output (filtered_feature_bc_matrix) to finished Seurat/Scanpy objects for analysis.

- Gene_TCR_Analysis
Runs various analyses using GEX and TCR for both mouse and human. Uses primarily Scanpy objects. 

- AUC_ROC_MarkerAnalyses
Runs various analyses using GEX, TCR, and COMET data on both mouse and human data. Uses primarily Seurat
objects.

- MachineLearning
Machine Learning analysis for Figure 2 of the paper.
