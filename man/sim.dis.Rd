\name{sim.dis}
\alias{sim.dis}
\title{Similarity Statistic Subroutine}
\description{
    Internal subroutine called by \code{\link[MVR]{MeanVarReg}} subroutine.
    Computes similarity statistic with sampling variability under a standard Gaussian reference distibution.
    Offers parallel computation option for computational efficiency.
}
\usage{
    sim.dis(data, k, B, parallel, conf, verbose)
}
\arguments{
  \item{data}{Standardized \code{numeric} data \code{matrix}, where points to cluster are by rows (usually samples), 
              or an object that can be coerced to such a \code{matrix} 
              (such as a \code{vector} or a \code{data.frame} with all \code{numeric} entries).
              Missing values are not allowed. \code{NaN} or \code{Inf} are not allowed.}
  \item{k}{Positive \code{integer} scalar of the fixed number of clusters.}
  \item{B}{Positive \code{integer} scalar of the number of Monte Carlo replicates of the inner loop of the gap statistic function.}
  \item{parallel}{\code{logical} scalar. Is parallel computing to be performed? Optional, defaults to \code{FALSE}.}           
  \item{conf}{\code{list} of parameters for cluster configuration, 
              passed from \code{\link[MVR]{MeanVarReg}} subroutine for parallel computing.
              Inputs for R package \pkg{parallel} function \code{makeCluster} (R package \pkg{parallel}) for cluster setup.
              See details for usage.}
  \item{verbose}{\code{logical} scalar. Is the output to be verbose?}
}
\details{
    See details in \code{\link[MVR]{mvr}} or \code{\link[MVR]{mvrt.test}}.
}
\value{
    \item{gapk}{\code{numeric} scalar of the similarity statistic value.}
    \item{sk}{\code{numeric} scalar of the standard error of the similarity statistic value.}
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
    Internal subroutine of \code{\link[MVR]{MeanVarReg}}, not to be called by end-user.
}
\seealso{
    \itemize{
        \item\code{\link[MVR]{MeanVarReg}} Mean-Variance Regularization Core Function.
        \item\code{\link[MVR]{withinsumsq}} Within-Cluster Sum of Squares Distances.
    }
}
\keyword{Gap Statistic}
\keyword{Reference Distribution}
\keyword{Monte-Carlo Replication}
\keyword{Parallel Programming}
