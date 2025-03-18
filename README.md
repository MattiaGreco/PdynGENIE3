# PdynGENIE3
This repository contains an implementation of the dynGENIE3 algorithm for network inference from time series data, with additional features.
The original algorithm was developed by "Huynh-Thu V., Geurts P., dynGENIE3: dynamical GENIE3 for the inference of gene networks from time series expression data. Sci Rep 8, 3384 (2018).". 
The following extension was created by Matia Greco and Olivier C. Martin. The new features are explained in full detail in "Greco M., Ricci-Tersenghi F., Martin O.C., Boosting reliability when inferring interactions from time series data in gene regulatory networks ".

In the following repository one can find:

-tutorial_PdynGENIE3.R : Example use of the code with priors;

-UsefulFuncitons.R : A few necessary functions to compute correlation matrices, decay rates and to do I/O;

-PdynGENIE3.c : Modified C code with the implementation of new features, in particular in this implementation the random forest construction is biased using priors and additional quantities are computed. In particular the new code returns: 
    -n_eff: the effective number of effective drivers of a given target;
    -alphas: the new estimates of the decay rates;
    -forest_quality: Average square error on the test set of the random forest regressor (1 number per gene);
    -square_error: Average square error on the test set of a give tree normalized by the variance at the root of the tree;
    -ls_size: size of the learning set used in the tree construction;

-PerformanceEstimation : Codes used to compute AUROC and AUPRC on the data;

-TestData_InOut1 : 5 graphs of 100 nodes with in-out degree of 1 returnd from Gene Net Weaver.



