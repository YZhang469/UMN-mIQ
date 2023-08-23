# R Code for "Modified Interactive Q-learning for Attenuating the Impact of Model Misspecification with Treatment Effect Heterogeneity"
Authors: Yuan Zhang, David M. Vock, Megan E. Patrick, Thomas A. Murray

Code is written by Y. Zhang. Please address all questions or comments to yuan.zhang@pennmedicine.upenn.edu.

This repository contains R code for the simulations and data analyses conducted in the manuscript. The functions for interactive Q-learning are modified based on the R package `iqLearn` developed by Linn et al. (2015).

## Simulation study

The document [simulation.R](https://github.com/YZhang469/UMN-mIQ/blob/main/simulation/simulation.R) contains the code for the preliminary simulation study (Supplementary Materials Appendix C) and the main simulation study (Section 5). In the simulation, we use a test dataset with all counterfactual outcomes to assess the performance of each model. Therefore, [helpFunctions.R](https://github.com/YZhang469/UMN-mIQ/blob/main/simulation/helpFunctions.R) and [functions.R](https://github.com/YZhang469/UMN-mIQ/blob/main/simulation/functions.R) are required.

## Application
The document [analysis.R](https://github.com/YZhang469/UMN-mIQ/blob/main/application/analysis.R) contains the code for the data analysis (results shown in Section 6).

## Reference
* Linn, K. A., Laber, E. B., & Stefanski, L. A. (2015). iqLearn: Interactive Q-Learning in R. Journal of Statistical Software, 64(1), 1–25. https://doi.org/10.18637/jss.v064.i01
