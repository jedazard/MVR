=======
MVR
=======
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

MVR is Open Source / Free Software, and is freely available under the GNU General Public License, version 3.

==========
References
==========
The companion papers can be accessed here:

Comput. Statist. Data Anal. (2012):
http://www.sciencedirect.com/science/article/pii/S0167947312000321

ASA-IMS JSM Proceedings (2011): 
https://www.amstat.org/membersonly/proceedings/2011/papers/302266_68145.pdf

ASA-IMS JSM Proceedings (2010): 
https://www.amstat.org/membersonly/proceedings/2010/papers/309104_62376.pdf

See also below on how to cite the package with the R command: citation("MVR").

=========================
Documentation and Manual: 
=========================
All the codes are in the R folder and a manual (MVR.pdf) details the end-user (and internal) functions. At this stage and for simplicity, there are only 2 end-user function, 4 end-user diagnostic and plotting functions and 2 end-user datasets (synthetic and real). See the "MVR-package" introduction section of the manual for more details and examples of use.

=============
Installation: 
=============
MVR 1.20.0 was built under R version 3.1.2 (2014-10-31).
Installation has been tested on Windows, Linux and Mac platforms.
To install the current software, download the version from the CRAN repository:

http://cran.r-project.org/web/packages/MVR/index.html


To install the development version 1.21.0, simply type:

library(devtools)

devtools::install_github("jedazard/MVR")

=============
Usage: 
=============
To load the MVR library in an R session and start using it:

library("MVR")

MVR.news()

citation("MVR")

etc...
