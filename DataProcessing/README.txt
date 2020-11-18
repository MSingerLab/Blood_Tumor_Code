This directory, DataProcessing, holds a few scripts to process the single-cell RNA sequencing data.
These scripts should be used after running Cell Ranger. The order is as follows:

1.) Run the 'clone_pipeline.py' script from the command line using the VDJ output from Cell Ranger as input.
This script takes in two samples (one from blood and one from tumor) as input.
This produces two csv files that contain the TCR data, and will be used in the next script.
Note that this script is run using Python 2 and not Python 3.

2.) Run the 'filtered_feature_bc_matrix' through '10X_to_Filtered_rds.R' script.
This script, using the gene lists in the GeneLists directory, filters and clusters all the cells in the
scRNA data on a sample-by-sample basis.

3.) Run the 'Define_Matching.R' script to determine matching status and clone size for each T cell,
using the filtered data objects (the outputs of Step 2).

4.) Integrate the data (IntegrateMouseData.R, IntegrateHumanData.R, and Merge_K409_LN_Tumor.R) to
produce integrated Seurat/rds objects. 

5.) The filtered objects can directly be used in the AUC_ROC_MarkerAnalyses directory.

6.) To run Python analyses, run the filtered data objects through 'Convert_seurat_csv.R' to turn
the Seurat objects to csv files, and then run 'Convert_scanpy_csv.ipynb' to convert the csv files
into Scanpy objects.

7.) Run the Scanpy objects through the DataAnalysis directory and scripts. 
