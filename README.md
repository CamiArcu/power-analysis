# Power Analysis

The purpose of Power Analysis is to evaluate the ability of a given program/method to solve a task. In this particular case, the task is detecting differentially expressed genes in a single-cell analysis. To evaluate this task, the program uses a real data set (composed by the gene expression of a single cell subtype) to simulate a new group with expression changes in only some predifined genes (which we will categorize as positive). Then it executes a differential expression analysis between both groups, using a user specified method, which will correctly detect some positives genes, while not others, thus obtaining the number of true positives, true negatives, false positives and false negatives tests. With these approach, the program calculates the power (sensitivity), the alpha (1 - specificity) or the area under the curve (AUC) of the selected method.


As in single-cell analysis gene expression changes between groups could be detected at the average expression value or the proportion of cells in which it is expressed, the program simulatse both types of changes allowing genes to only change their average expression, proportion of cells expressing it or both, in any direction.
Power, alpha and AUC are calculated by the program for each of the differential expression change types and the change-magnitus could be user specified.

Also, as the quality of the single cell analysis is much influenced by the number of initial cells to be evaluated and its library size, the program calculates each of the metrics for multiple user specified total cells amount and average library size.


All these features could be visually integrated to determine the best method, sequenced cell number and expected library size to detect average expression changes, differences in propotion of expressing cells or both.

