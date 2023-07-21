# R Code for "Modified Interactive Q-learning for Attenuating the Impact of Model Misspecification with Treatment Effect Heterogeneity"
Authors: Yuan Zhang, David M. Vock, Megan E. Patrick, Thomas A. Murray

Code was written by Y. Zhang. Please address all questions or comments to yuan.zhang@pennmedicine.upenn.edu.

This repository contains R code for the simulations conducted in the manuscript.

The functions for interactive Q-learning are modified based on Linn et al. (2015). The `iqLearn` functions need to be modified because, in the simulation, we use a test dataset with all counterfactual outcomes to assess the performance of each model.

## Reference
* Linn, K. A., Laber, E. B., & Stefanski, L. A. (2015). iqLearn: Interactive Q-Learning in R. Journal of Statistical Software, 64(1), 1–25. https://doi.org/10.18637/jss.v064.i01
