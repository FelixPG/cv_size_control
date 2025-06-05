# cv_size_control
Repository for the codes and data associated with the Proulx-Giraldeau et al., 2025 paper. 


- Movies.
  Folder containing the data analyzed in the paper. This folder contains .csv files outputed by the Cell-ACDC segmentation and tracking software (Padovani, 2022) created from the analysis of phase-contrast microscopy movies of live cells. This include 2 replicate experiments of cln3 deletion cells, 2 replicate experiments for cln3whi5 double deletion mutant cells, and 2 replicate experiments for WT cells grown in glucose (data provided by Yagya Chadha and Kurt Schmoller). Additionally, this folder contains the data provided by Yuping Chen and Bruce Futcher from their Mol. Cell 2020 paper as presented in Figure 1 of our publication. Lastly, this also contains the data for the asymmetry analysis done in Figure 5 of our publication.

- main_downsteam_movie_analysis.ipynb 
  Main Python notebook for analyzing the cln3 and cln3whi5 mutants described in this publication. Data is contained in the Movies folder.

- downstream_movie_analysis_WT.ipynb
  Variant of the Python notebook for analyzing movie data. This one pertains to the WT data in glucose from the Schmoller lab which is included in the Movies folder.

- compare_mutants.ipynb
  Python notebook for analyzing the movie data and compare mutants against each other.

- chen_futcher_data_processing.ipynb
  Python notebook for processing and analyzing the data from the Chen, 2020 paper included in the Movies folder, and producing the plots of Figure 1 of this publication. 

- cell_cycle_model.m 
  Matlab script that simulates the model from Chandler-Brown, 2017. Adapted for the purposes of this publication. This script outputs a Matlab workspace that can be further analyzed by other scripts.
  
- analyze_run.mlx 
  Matlab notebook to analyze a simulation from cell_cycle_model.m and produce various plots for the figures 2, 3 and 4 of this publication.

- compare_mutants.mlx
  Matlab notebook to analyze many simulations from cell_cycle_model.m and produce plots from the Figures 4 and 5 of this publication. 
  
- parameter_sensitivity.m
  Matlab script that computes the numerical Hessian of the mean, the variance and the CV of a distribution of your choice from cell_cycle_model.m, and outputs the eigenvalues and eigenvectors of the matrix, the first one being the gradient of the measure under study.

- brewermap.m
  Matlab function for plotting.

- shadederrormap.m
  Matlab function for plotting. 
