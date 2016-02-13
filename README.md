### General Remarks

CRAN downloads  October 1, 2012, 
the month the [RStudio CRAN mirror](http://cran-logs.rstudio.com/) 
started publishing logs:
[![](http://cranlogs.r-pkg.org/badges/grand-total/MVR)](http://cran.rstudio.com/web/packages/MVR/index.html)

CRAN downloads in the last month:
[![](http://cranlogs.r-pkg.org/badges/last-month/MVR)](http://cran.rstudio.com/web/packages/MVR/index.html)

CRAN downloads in the last week:
[![](http://cranlogs.r-pkg.org/badges/last-week/MVR)](http://cran.rstudio.com/web/packages/MVR/index.html)


===============
### Description

MVR (Mean-Variance Regularization) is a non-parametric method for joint adaptive mean-variance regularization 
and variance stabilization of high-dimensional data. It is suited for handling difficult problems posed by 
high-dimensional multivariate datasets (_p_ >> _n_ paradigm), such as in omics-type data, 
among which are that the variance is often a function of the mean, variable-specific estimators of variances are not reliable, 
and tests statistics have low powers due to a lack of degrees of freedom.

Key features include:

1. Normalization and/or variance stabilization of the data

2. Computation of mean-variance-regularized _t_-statistics (_F_-statistics to come)

3. Generation of diverse diagnostic plots

4. Computationally efficient implementation using C/C++ interfacing and an option for parallel
computing to enjoy a fast and easy experience in the R environment

See also below the package news with the R command: `MVR.news()`.


============================
### Documentation and Manual

All the codes are in the R folder and a manual (MVR.pdf) details the end-user (and internal) functions. 
At this stage and for simplicity, there are only 2 end-user function, 4 end-user diagnostic 
and plotting functions and 2 end-user datasets (synthetic and real). 
See the "MVR-package" introduction section of the manual for more details and examples.


==============
### References

Open access to companion papers (accepted for publication):

- [Comput. Statist. Data Anal. (2012)](http://www.sciencedirect.com/science/article/pii/S0167947312000321).
The Official Journal of the International Association for Statistical Computing.

- [JSM (2011)](https://www.amstat.org/membersonly/proceedings/2011/papers/302266_68145.pdf). 
The ASA Proceedings of the annual Joint Statistical Meetings (Miami, FL, USA).

- [JSM (2010)](https://www.amstat.org/membersonly/proceedings/2010/papers/309104_62376.pdf). 
The ASA Proceedings of the annual Joint Statistical Meetings (Vancouver, BC, Canada).


===========
### License

MVR is open source / free software, licensed under the GNU General Public License, version 3 (GPLv3), 
sponsored by the [Free Software Foundation](http://www.fsf.org/). To view a copy of this license, visit 
[GNU Free Documentation License](http://www.gnu.org/licenses/gpl-3.0.html).


================
### Installation

* To install MVR from CRAN, simply download and install the current version (1.30.3) from the CRAN repository:

```{r}
install.packages("MVR")
```

* Alternatively, you can install the most up-to-date development version (1.30.3) from GitHub, using devtools:

```{r}
install.packages("devtools")
library("devtools")
devtools::install_github("jedazard/MVR")
```

================
### Requirements

MVR 1.30.3 requires R-3.0.2 (2013-09-25). It was built and tested under R-devel (2015-12-22 r69809) and Travis CI. 

Installation has been tested on Windows, Linux, OSX and Solaris platforms. 

See Travis CI build result:

[![Build Status](https://travis-ci.org/jedazard/MVR.png?branch=master)](https://travis-ci.org/jedazard/MVR)

See CRAN checks:

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MVR)](https://cran.r-project.org/web/checks/check_results_MVR.html).


=========
### Usage

* To load the MVR library in an R session and start using it:

```{r}
library("MVR")
```

* Check the package news with the R command:

```{r}
MVR.news()
```

* Check on how to cite the package with the R command:

```{r}
citation("MVR")
```

etc...
