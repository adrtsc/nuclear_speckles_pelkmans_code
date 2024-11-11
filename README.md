# nuclear_speckles_pelkmans_code
Code resources associated with the publication "Phosphorylation of a nuclear condensate regulates cohesion and mRNA retention "

This repository contains two folders:

## transcript_feature_calculation:
Scripts used to generate transcript features (folder transcript_feature_calculation). The following files have been deposited:
  - get_transcript_features_20240724.py: Computes transcript features given annotations and a FASTA record dict
  - get_transcript_sequences_20231024.py: Collects transcript sequences from FASTA record given transcript annotations
  - run_transcript_features: Script to run get_transcript_features_20240724.py on a slurm cluster

## MLR_scripts
Scripts used to run multilinear regression model to predict nulcear speckle enrichment based on transcript features. Scripts have a pre-fix indicating the order in which they are run:
  
  - 00_prepare_data.py: Script to collect and prepare data before running MLR models.
  - 00_run_prepare_data.py: Script to run step 00 on slurm cluster.
  - 01_MLR_10_fold_cv.py: Script to run MLR with 10-fold cross validation. 
  - 01_run_MLR.py: Script to run step 01 on slurm cluster.
  - 01b_MLR_10_fold_cv_Huber.py: Script to run MLR with 10-fold cross validation using a Huber linear regressor. 
  - 01b_run_MLR_Huber.py: Script to run step 01b on slurm cluster.
  - 02_MLR_combinatorial.py: Script to run MLR model on different subsets of feature values.
  - 02_run_MLR_combinatorial.py: Script to run step 02 on slurm cluster.
  - 03_MLR_top_loadings.py: Script to analyse top PC loadings in PCs with highest absolute coefficients in MLR.
  - 03_run_MLR_top_loadings.py: Script to run step 03 on slurm cluster.
