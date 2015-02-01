\name{MeanVarReg}
\alias{MeanVarReg}
\title{Mean-Variance Regularization Core Subroutine}
\description{
    Core subroutine of \code{\link[MVR]{mvr}} and \code{\link[MVR]{mvrt.test}} functions for
    mean-variance regularization and variance stabilization, and the computation of regularized test-statistics. 
    Returns optimal cluster configuration. 
    Internal function called by \code{\link[MVR]{mvr}} and internally by \code{\link[MVR]{mvrt}}.
}
\usage{
    MeanVarReg(data, nc.min, nc.max, probs, B, parallel, conf, verbose)
}
\arguments{
  \item{data}{\code{numeric} data \code{matrix}, where points to cluster are by rows (usually samples), 
              or an object that can be coerced to such a \code{matrix} 
              (such as a \code{vector} or a \code{data.frame} with all \code{numeric} entries).
              Missing values (\code{NA}), NotANumber values (\code{NaN}),
              or Infinite values (\code{Inf}) are not allowed.}
  \item{nc.min}{Positive \code{integer} scalar of the minimum number of clusters.}
  \item{nc.max}{Positive \code{integer} scalar of the maximum number of clusters.}
  \item{B}{Positive \code{integer} scalar of the number of Monte Carlo replicates 
           of the inner loop of the sim statistic function.}
  \item{probs}{\code{numeric} \code{vector} of probabilities for quantile diagnostic plots.}
  \item{parallel}{\code{logical} scalar. Is parallel computing to be performed? Optional, defaults to \code{FALSE}.}           
  \item{conf}{\code{list} of parameters for cluster configuration, passed from \code{\link[MVR]{mvr}} 
              or \code{\link[MVR]{mvrt.test}} functions for parallel computing.
              Inputs for R package \pkg{parallel} function \code{makeCluster} (R package \pkg{parallel}) for cluster setup.
              See details for usage.}
  \item{verbose}{\code{logical} scalar. Is the output to be verbose?}
}
\details{
        See details in \code{\link[MVR]{mvr}} and \code{\link[MVR]{mvrt.test}}.
}
\value{
    \item{membership}{\code{numeric} \code{vector} of cluster membership of each variable}
    \item{nc}{Positive \code{integer} scalar of number of clusters found in optimal cluster configuration}
    \item{gap}{\code{numeric} \code{vector} of the similarity statistic values}
    \item{sde}{\code{numeric} \code{vector} of the standard errors of the similarity statistic values}
    \item{mu.quant}{\code{numeric} \code{matrix} (\code{nc.max} - \code{nc.min} + 1) x (length(\code{probs})) of quantiles of means}
    \item{sd.quant}{\code{numeric} \code{matrix} (\code{nc.max} - \code{nc.min} + 1) x (length(\code{probs})) of quantiles of standard deviations}
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
        \item "Jean-Eudes Dazard, Ph.D." \email{jean-eudes.dazard@case.edu}
        \item "Hua Xu, Ph.D." \email{huaxu77@gmail.com}
        \item "Alberto Santana, MBA." \email{ahs4@case.edu}
    }
    Maintainer: "Jean-Eudes Dazard, Ph.D." \email{jean-eudes.dazard@case.edu}
    
    Acknowledgments: This project was partially funded by the 
    National Institutes of Health (P30-CA043703 to J-E.DAZARD).
}
\note{
    Internal function not to be called by end-user.
}
\seealso{
    \itemize{
        \item\code{\link[MVR]{mvr}} Mean-Variance Regularization and Variance Stabilization.
        \item\code{\link[MVR]{mvrt.test}} Mean-Variance Regularized t-Test Statistic with Significance.
        \item\code{\link[MVR]{mvrt}} Mean-Variance Regularized t-test statistic.
    }
}
\keyword{mean-variance estimators}
\keyword{regularization}
\keyword{variance stabilization}
\keyword{normalization}
\keyword{regularized test statistics}
\keyword{Parallel Programming}
