# Power Analysis

The program is an adaptation of the _POWSC_ algorithm implemented by [Suke](https://github.com/suke18/POWSC)

The purpose of power-analysis is to evaluate the ability of a given program/method to solve a task. In this particular case, the task is detecting differentially expressed genes in a single-cell analysis. To evaluate this task, the program uses a real data set (composed by the gene expression of a single cell subtype) to simulate a new group with expression changes in only some predifined genes (which we will categorize as positive). Then it executes a differential expression analysis between both groups, using a user specified method, which will correctly detect some positives genes, while not others, thus obtaining the number of true positives, true negatives, false positives and false negatives tests. With these approach, the program calculates the power (sensitivity), the alpha (1 - specificity) or the area under the curve (AUC) of the selected method.


As in single-cell analysis gene expression changes between groups could be detected at the average expression value or the proportion of cells in which it is expressed, the program simulatse both types of changes allowing genes to only change their average expression, proportion of cells expressing it or both, in any direction.
Power, alpha and AUC are calculated by the program for each of the differential expression change types and the change-magnitus could be user specified.

Also, as the quality of the single cell analysis is much influenced by the number of initial cells to be evaluated and its library size, the program calculates each of the metrics for multiple user specified total cells amount and average library size.


All these features could be visually integrated to determine the best method, sequenced cell number and expected library size to detect average expression changes, differences in propotion of expressing cells or both.

## How to run It

Dependencies needed:\n
POWSC. Requires R > =  4.1
tidyverse, purrr, viridis, scran, DESeq2, zinbwave, BiocParallel, furrr, RhpcBLASctl, gridExtra, pROC, limma

The power-analysis was thought to be run from terminal, executing the **PowerAnalysis_genefilter_DE1_DE2_comb.sh** file.\n 

This file will execute the **PowerAnalysis_genefilter_dpig_mu_comb.R** script with the setted:

- _j_ vector of total cells (half of total cells for each evaluated group),
- _mu_ vector of absolut log2FC of average expression between evaluated groups,
- _k_ vector of absolute difference of proportion of cells expressing the gene between evaluated groups and
- _i_ vector of run number, wich length will correspond t the times the analysis will be iterated.

Similar conditions are displayed on the **PowerAnalysis_genefilter_DE1_DE2_comb_diffUMICell.sh** which was tought to evaluate different number of total UMIs on the dataset.

All other perameters remains the same but the program will incoparte:

- _l_ vector of tranformed library sizes


User can change this parameters changinh the bash files.

All the analysis will be executed evaluating the DESeq2-zinbwave method of differential expression analysis, but user can change this looking at the **PowerAnalysis_genefilter_dpig_mu_comb.R** or **PowerAnalysis_genefilter_dpig_mu_comb_differentUMIcell.R** files, changing the _DE_Method_ parameter for all other possibles:

- MAST,
- SC2P,
- DESeq2, 
- ScDESeq2,  
- SeuratDESeq2, 
- DESeq2zw: DESeq2-zinbwave


The "pkill R" command and the structure of the bash files, which divies them by blocks was implemented to avoid memory issues on my personall computer but could be chages on local repositories according the user devices capabilityies.

## Analysis of te results

Examples of the analysis of the results are stored at the **Analizing_PowerAnalysis_genefilter_dpig_mu_comb_220125.R** and **Analizing_PowerAnalysis_genefilter_dpig_mu_comb_220225.R** files. It is possible to locally modify these files to obtain heatmaps, scatter pplots and AUC curves of the evaluated parameters.



 
