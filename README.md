=======
MVR
=======
Mean-Variance Regularization (MVR) is a non-parametric method for joint adaptive mean-variance regularization and variance stabilization of high-dimensional data.
It is suited for handling difficult problems posed by high-dimensional multivariate datasets (p  n
paradigm), such as in omics-type data, among which are that the variance is often a function of the
mean, variable-specific estimators of variances are not reliable, and tests statistics have low powers
due to a lack of degrees of freedom.
Key features include:
1. Normalization and/or variance stabilization of the data
2. Computation of mean-variance-regularized t-statistics (F-statistics to come)
3. Generation of diverse diagnostic plots
4. Computationally efficient implementation using C/C++ interfacing and an option for parallel
computing to enjoy a fast and easy experience in the R environment

=========================
Documentation and Manual: 
=========================
All the codes are in the R folder and a manual (MVR.pdf) details the end-user (and internal) functions. At this stage and for simplicity, there are only two end-user function, 4 end-user diagnostic and plotting functions and 2 end-user datasets (synthetic and real). See the "MVR-package" introduction section of the manual for more details and examples of use.

=============
Installation: 
=============
The latest R version 3.1.2 (2014-10-31) is required.
Installation has been tested on Windows, Linux and Mac platforms.
To install the software and load the MVR library in an R session, simply type:

library(devtools)

devtools::install_github("jedazard/MVR")

library("MVR")
