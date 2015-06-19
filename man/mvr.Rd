\name{mvr}

\alias{mvr}

\title{Function for Mean-Variance Regularization and Variance Stabilization}

\description{
    End-user function for Mean-Variance Regularization (MVR) and Variance Stabilization by similarity statistic
    under sample group homoscedasticity or heteroscedasticity assumptions.

    Return an object of class "\code{mvr}". Offers the option of parallel computation for improved efficiency.
}

\usage{
    mvr(data,
        block = rep(1,nrow(data)),
        tolog = FALSE,
        nc.min = 1,
        nc.max = 30,
        probs = seq(0, 1, 0.01),
        B = 100,
        parallel = FALSE,
        conf = NULL,
        verbose = TRUE)
}

\arguments{
  \item{data}{\code{numeric} \code{matrix} of untransformed (raw) data,
              where samples are by rows and variables (to be clustered) are by columns,
              or an object that can be coerced to such a \code{matrix} (such as a \code{numeric} \code{vector}
              or a \code{data.frame} with all \code{numeric} columns).
              Missing values (\code{NA}), NotANumber values (\code{NaN}) or Infinite values (\code{Inf}) are not allowed.}
  \item{block}{\code{character} or \code{numeric} \code{vector} or \code{factor} grouping/blocking variable of length the sample size.
               Defaults to single group situation (see details).}
  \item{tolog}{\code{logical} scalar. Is the data to be log2-transformed first? Optional, defaults to \code{FALSE}.
             Note that negative or null values will be changed to 1 before taking log2-transformation.}
  \item{nc.min}{Positive \code{integer} scalar of the minimum number of clusters, defaults to 1}
  \item{nc.max}{Positive \code{integer} scalar of the maximum number of clusters, defaults to 30}
  \item{probs}{\code{numeric} \code{vector} of probabilities for quantile diagnostic plots. Defaults to \code{seq}(0, 1, 0.01).}
  \item{B}{Positive \code{integer} scalar of the number of Monte Carlo replicates of the inner loop
           of the sim statistic function (see details).}
  \item{parallel}{\code{logical} scalar. Is parallel computing to be performed? Optional, defaults to \code{FALSE}.}
  \item{conf}{\code{list} of parameters for cluster configuration.
              Inputs for R package \pkg{parallel} function \code{makeCluster} (R package \pkg{parallel}) for cluster setup.
              Optional, defaults to \code{NULL}. See details for usage.}
  \item{verbose}{\code{logical} scalar. Is the output to be verbose? Optional, defaults to \code{TRUE}.}
}

\details{
    Argument \code{block} is a \code{vector} or a \code{factor} grouping/blocking variable. It must be of length sample size
    with as many different \code{character} or \code{numeric} values as the number of levels or sample groups.
    It defaults to single group situation, i.e. under the assumption of equal variance between sample groups.
    All group sample sizes must be greater than 1, otherwise the program will stop.

    Note that argument \code{B} is internally reset to \code{conf$cpus}*\code{ceiling}(\code{B}/\code{conf$cpus}) in case the parallelization
    is used (i.e. \code{conf} is non \code{NULL}), where \code{conf$cpus} denotes the total number of CPUs to be used (see below).

    Argument \code{nc.max} currently defaults to 30. Empirically, we found that this is enough for most datasets tested.
    This depends on (i) the dimensionality/sample size ratio \eqn{\frac{p}{n}}, (ii) the signal/noise ratio, and
    (iii) whether a pre-transformation has been applied (see \cite{Dazard, J-E. and J. S. Rao (2012)} for more details).
    See the cluster diagnostic function \code{\link[MVR]{cluster.diagnostic}} for more details, whether larger values of \code{nc.max} may be required.

    To run a parallel session (and parallel RNG) of the MVR procedures (\code{parallel}=\code{TRUE}), argument \code{conf}
    is to be specified (i.e. non \code{NULL}). It must list the specifications of the folowing parameters for cluster configuration:
    "names", "cpus", "type", "homo", "verbose", "outfile". These match the arguments described in function \code{makeCluster}
    of the R package \pkg{parallel}. All fields are required to properly configure the cluster, except for "names" and "cpus",
    which are the values used alternatively in the case of a cluster of type "SOCK" (socket), or in the case of a cluster
    of type other than "SOCK" (socket), respectively.
    \itemize{
        \item "names": \code{names} : \code{character} vector specifying the host names on which to run the job.
              Could default to a unique local machine, in which case, one may use the unique host name "localhost".
              Each host name can potentially be repeated to the number of CPU cores available on the corresponding machine.
        \item "cpus": \code{spec} : \code{integer} scalar specifying the total number of CPU cores to be used
              across the network of available nodes, counting the workernodes and masternode.
        \item "type": \code{type} : \code{character} vector specifying the cluster type ("SOCK", "PVM", "MPI").
        \item "homo": \code{homogeneous} : \code{logical} scalar to be set to \code{FALSE} for inhomogeneous clusters.
        \item "verbose": \code{verbose} : \code{logical} scalar to be set to \code{FALSE} for quiet mode.
        \item "outfile": \code{outfile} : \code{character} vector of the output log file name for the workernodes.
    }
    
    The actual creation of the cluster, its initialization, and closing are all done internally.
    In addition, when random number generation is needed, the creation of separate streams of parallel RNG per node
    is done internally by distributing the stream states to the nodes (For more details see function \code{makeCluster}
    (R package \pkg{parallel}) and/or \url{http://www.stat.uiowa.edu/~luke/R/cluster/cluster.html}.
}

\value{
    \item{Xraw}{\code{numeric} \code{matrix} of original data.}
    \item{Xmvr}{\code{numeric} \code{matrix} of MVR-transformed data.}
    \item{centering}{\code{numeric} \code{vector} of centering values
                     for standardization (cluster mean of pooled sample mean).}
    \item{scaling}{\code{numeric} \code{vector} of scaling values
                   for standardization (cluster mean of pooled sample std dev).}
    \item{MVR}{\code{list} (of size the number of groups) containing for each group:
        \itemize{
            \item membership {\code{numeric} \code{vector} of cluster membership of each variable}
            \item nc {Positive \code{integer} scalar of number of clusters found in optimal cluster configuration}
            \item gap {\code{numeric} \code{vector} of the similarity statistic values}
            \item sde {\code{numeric} \code{vector} of the standard errors of the similarity statistic values}
            \item mu.std {\code{numeric} \code{matrix} (\code{K} x p) of the vector of standardized means by groups (rows), 
                          where \code{K} = \#groups and \code{p} = \#variables}
            \item sd.std {\code{numeric} \code{matrix} (\code{K} x p) of the vector of standardized standard deviations by groups (rows),
                          where \code{K} = \#groups and \code{p} = \#variables}
            \item mu.quant {\code{numeric} \code{matrix} (\code{nc.max} - \code{nc.min} + 1) x (length(\code{probs})) of quantiles of means}
            \item sd.quant {\code{numeric} \code{matrix} (\code{nc.max} - \code{nc.min} + 1) x (length(\code{probs})) of quantiles of standard deviations}
        }
    }
    \item{block}{Value of argument \code{block}.}
    \item{tolog}{Value of argument \code{tolog}.}
    \item{nc.min}{Value of argument \code{nc.min}.}
    \item{nc.max}{Value of argument \code{nc.max}.}
    \item{probs}{Value of argument \code{probs}.}
}

\references{
    \itemize{
        \item Dazard J-E., Hua Xu and J. S. Rao (2011). "\emph{R package MVR for Joint Adaptive Mean-Variance Regularization and Variance Stabilization.}"
              In JSM Proceedings, Section for Statistical Programmers and Analysts. Miami Beach, FL, USA: American Statistical Association IMS - JSM, 3849-3863.
        \item Dazard J-E. and J. S. Rao (2012). "\emph{Joint Adaptive Mean-Variance Regularization and Variance Stabilization of High Dimensional Data.}"
              Comput. Statist. Data Anal. 56(7):2317-2333.
    }
}

\author{
    \itemize{
        \item "Jean-Eudes Dazard, Ph.D." \email{jxd101@case.edu}
        \item "Hua Xu, Ph.D." \email{huaxu77@gmail.com}
        \item "Alberto Santana, MBA." \email{ahs4@case.edu}
    }
    Maintainer: "Jean-Eudes Dazard, Ph.D." \email{jxd101@case.edu}

    Acknowledgments: This project was partially funded by the
    National Institutes of Health (P30-CA043703 to J-E.DAZARD).
}

\note{
    End-user function.
}

\seealso{
    \itemize{
        \item \code{makeCluster} (R package \pkg{parallel}).
        \item \code{justvsn} (R package \pkg{vsn}) Variance stabilization and calibration for microarray data \cite{Huber, 2002}
    }
}

\examples{
#===================================================
# Loading the library and its dependencies
#===================================================
library("MVR")

\dontrun{
    #===================================================
    # MVR package news
    #===================================================
    MVR.news()

    #================================================
    # MVR package citation
    #================================================
    citation("MVR")

    #===================================================
    # Loading of the Synthetic and Real datasets
    # (see description of datasets)
    #===================================================
    data("Synthetic", "Real", package="MVR")
    ?Synthetic
    ?Real
}

#===================================================
# Mean-Variance Regularization (Synthetic dataset)
# Single-Group Assumption
# Assuming equal variance between groups
# Without cluster usage
#===================================================
nc.min <- 1
nc.max <- 10
probs <- seq(0, 1, 0.01)
n <- 10
mvr.obj <- mvr(data = Synthetic,
               block = rep(1,n),
               tolog = FALSE,
               nc.min = nc.min,
               nc.max = nc.max,
               probs = probs,
               B = 100,
               parallel = FALSE,
               conf = NULL,
               verbose = TRUE)

\dontrun{
    #===================================================
    # Examples of parallelization below with 
    # a SOCKET or MPI cluster configuration
    #===================================================
    # 1- WINDOWS multicores PC with SOCKET communication
    #    With a 2-Quad (8-CPUs) PC
    #===================================================
    if (.Platform$OS.type == "windows") {
        cpus <- detectCores()
        conf <- list("names" = rep("localhost", cpus),
                     "cpus" = cpus,
                     "type" = "SOCK",
                     "homo" = TRUE,
                     "verbose" = TRUE,
                     "outfile" = "")
    }
    #===================================================
    # 2- LINUX multinodes cluster with SOCKET communication
    #    with 4-nodes (32-CPUs) cluster
    #    with 1 masternode and 3 workernodes
    #    All hosts run identical setups
    #    Same number of core CPUs (8) per node
    #===================================================
    if (.Platform$OS.type == "unix") {
        masterhost <- Sys.getenv("HOSTNAME")
        slavehosts <- c("compute-0-0", "compute-0-1", "compute-0-2")
        nodes <- length(slavehosts) + 1
        cpus <- 8
        conf <- list("names" = c(rep(masterhost, cpus),
                                 rep(slavehosts, cpus)),
                     "cpus" = nodes * cpus,
                     "type" = "SOCK",
                     "homo" = TRUE,
                     "verbose" = TRUE,
                     "outfile" = "")
    }
    #===================================================
    # 3- LINUX multinodes cluster with MPI communication
    #    Here, a file named ".nodes" (e.g. in the home directory)
    #    must contain the list of nodes of the cluster
    #===================================================
    if (.Platform$OS.type == "unix") {
        hosts <- scan(file=paste(Sys.getenv("HOME"), "/.nodes", sep=""), 
                      what="", 
                      sep="\n")
        hostnames <- unique(hosts)
        nodes <- length(hostnames)
        cpus <-  length(hosts)/length(hostnames)
        conf <- list("cpus" = nodes * cpus,
                     "type" = "MPI",
                     "homo" = TRUE,
                     "verbose" = TRUE,
                     "outfile" = "")
    }
    #===================================================
    # Run:
    # Mean-Variance Regularization (Real dataset)
    # Multi-Group Assumption
    # Assuming unequal variance between groups
    #===================================================
    nc.min <- 1
    nc.max <- 30
    probs <- seq(0, 1, 0.01)
    n <- 6
    GF <- factor(gl(n = 2, k = n/2, len = n),
                 ordered = FALSE,
                 labels = c("M", "S"))
    mvr.obj <- mvr(data = Real,
                   block = GF,
                   tolog = FALSE,
                   nc.min = nc.min,
                   nc.max = nc.max,
                   probs = probs,
                   B = 100,
                   parallel = TRUE,
                   conf = conf,
                   verbose = TRUE)
    }
}

\keyword{Mean-Variance Estimators}
\keyword{Regularization}
\keyword{Variance Stabilization}
\keyword{Normalization}
\keyword{Parallel Programming}
\keyword{High Performance Computing}

