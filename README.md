### License

PRIMsrc is open source / free software, licensed under the GNU General Public License version 3 (GPLv3), 
sponsored by the [Free Software Foundation](http://www.fsf.org/). To view a copy of this license, visit 
[GNU Free Documentation License](http://www.gnu.org/licenses/gpl-3.0.html).


=============
### Downloads

CRAN downloads since October 1, 2012, 
the month the [RStudio CRAN mirror](http://cran-logs.rstudio.com/) 
started publishing logs:
[![](https://cranlogs.r-pkg.org/badges/grand-total/MVR)](https://CRAN.R-project.org/package=MVR)

CRAN downloads in the last month:
[![](https://cranlogs.r-pkg.org/badges/last-month/MVR)](https://CRAN.R-project.org/package=MVR)

CRAN downloads in the last week:
[![](https://cranlogs.r-pkg.org/badges/last-week/MVR)](https://CRAN.R-project.org/package=MVR)


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


================
### Installation

* To install [`MVR` from CRAN repository](https://CRAN.R-project.org/package=MVR), simply download and install the current version (1.31.0) from the CRAN repository:

```{r}
install.packages("MVR")
```

* Alternatively, you can install the most up-to-date development version (>= 1.32.0) of [`MVR` from GitHub repository](https://github.com/jedazard/MVR) using devtools, simply run:

```{r}
install.packages("devtools")
library("devtools")
devtools::install_github("jedazard/MVR")
```

================
### Requirements

MVR 1.31.0 requires R-3.0.2 (2013-09-25). It was built and tested under R-devel (2016-06-30 r70858) and Travis CI. 

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


===================
### Acknowledgments

Authors: 
   + Jean-Eudes Dazard, Ph.D. [(jean-eudes.dazard@case.edu)](jean-eudes.dazard@case.edu)
   + Hua Xu, Ph.D. [(huaxu77@gmail.com)](huaxu77@gmail.com)
   + Alberto Santana, MBA. [(ahs4@case.edu)](ahs4@case.edu)

Maintainers: 
   + Jean-Eudes Dazard, Ph.D. [(jean-eudes.dazard@case.edu)](jean-eudes.dazard@case.edu)

Funding/Provision/Help:   
   + This work made use of the High Performance Computing Resource in the Core Facility for Advanced Research Computing at Case Western Reserve University. 
   + This project was partially funded by the National Institutes of Health NIH - National Cancer Institute (P30-CA043703).


==============
### References
              
   + Dazard J-E. and J. S. Rao (2010). 
      *Regularized Variance Estimation and Variance Stabilization of High-Dimensional Data*. 
      In JSM Proceedings, Section for High-Dimensional Data Analysis and Variable Selection. 
      Vancouver, BC, Canada: American Statistical Association IMS - JSM, 5295-5309.
   + Dazard J-E., Hua Xu and J. S. Rao (2011). 
      *R package MVR for Joint Adaptive Mean-Variance Regularization and Variance Stabilization*. 
      In JSM Proceedings, Section for Statistical Programmers and Analysts. 
      Miami Beach, FL, USA: American Statistical Association IMS - JSM, 3849-3863.
   + Dazard J-E. and J. S. Rao (2012). 
      *Joint Adaptive Mean-Variance Regularization and Variance Stabilization of High Dimensional Data*. 
      Comput. Statist. Data Anal. 56(7):2317-2333.
