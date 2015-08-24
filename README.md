============
Description:
============
Mean-Variance Regularization (MVR) is a non-parametric method for joint adaptive mean-variance regularization and variance stabilization of high-dimensional data.
It is suited for handling difficult problems posed by high-dimensional multivariate datasets (p >> n
paradigm), such as in omics-type data, among which are that the variance is often a function of the
mean, variable-specific estimators of variances are not reliable, and tests statistics have low powers
due to a lack of degrees of freedom.
Key features include:

1. Normalization and/or variance stabilization of the data

2. Computation of mean-variance-regularized t-statistics (F-statistics to come)

3. Generation of diverse diagnostic plots

4. Computationally efficient implementation using C/C++ interfacing and an option for parallel
computing to enjoy a fast and easy experience in the R environment

See also below the package news with the R command: MVR.news().

========
License:
========
MVR is Open Source / Free Software, and is freely available under the GNU General Public License, version 3.

=========================
Documentation and Manual: 
=========================
All the codes are in the R folder and a manual (MVR.pdf) details the end-user (and internal) functions. At this stage and for simplicity, there are only 2 end-user function, 4 end-user diagnostic and plotting functions and 2 end-user datasets (synthetic and real). See the "MVR-package" introduction section of the manual for more details and examples.

==========
References
==========
CRAN release:
https://cran.r-project.org/web/packages/MVR/index.html


The companion papers can be accessed here:

- Comput. Statist. Data Anal. (2012):

http://www.sciencedirect.com/science/article/pii/S0167947312000321

- ASA-IMS JSM Proceedings (2011): 

https://www.amstat.org/membersonly/proceedings/2011/papers/302266_68145.pdf

- ASA-IMS JSM Proceedings (2010): 

https://www.amstat.org/membersonly/proceedings/2010/papers/309104_62376.pdf

=============
Requirements:
=============
MVR 1.30.3 requires R-3.0.2 (2013-09-25). It was built and tested under R-devel (2015-07-20 r68705) and Travis CI. 

Installation has been tested on Windows, Linux and OSX platforms. See for instance the 'CRAN Package Check Results' here:

https://cran.r-project.org/web/checks/check_results_MVR.html

=============
Installation: 
=============
- To install MVR from CRAN, download and install the current version (1.30.2) from the CRAN repository:

install.packages("MVR")

- Or, you can install the development version (1.30.3) from GitHub, using devtools:

install.packages(devtools)

library(devtools)

devtools::install_github("jedazard/MVR")

======
Usage: 
======
- To load the MVR library in an R session and start using it:

library("MVR")

- Check the package news with the R command:

MVR.news()

- Check on how to cite the package with the R command:

citation("MVR")

etc...
