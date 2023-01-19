# Source code of the paper: Explainable domain transfer of distant supervised cancer subtyping model via imaging-based rules extraction

FOLDERS:
1. DATA: folder of data (not included)
2. MatSurv-master: Code for log rank test and survival analysis 
   Cloned from https://github.com/aebergl/MatSurv
3. Subfunction (and S2GC.m script): Code for S2GC model
   Cloned from https://github.com/CLiu272/S2GC

Scripts:
a. transfer_clustering.m implements the analysis of the paper (DS-CS model)
b. best_parameters.m finds the best parameters for DS-CS model
c. concordance_index.m computes concordance index
d. findBestEnsTree.m fit the best ETM model
e. results-and-explainability.R run the explainability analysis and the rule extraction
f. kmean_clustering_comparison.m performs the naive k-mean clustering
