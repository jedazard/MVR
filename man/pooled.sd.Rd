\name{pooled.sd}
\alias{pooled.sd}
\title{Pooled Group Sample Standard Deviation Subroutine}
\description{
    Internal subroutine called by end-user function \code{\link[MVR]{mvr}} and internally by \code{\link[MVR]{MeanVarReg}} subroutine.
}
\usage{
    pooled.sd(x, block)
}
\arguments{
  \item{x}{\code{numeric} \code{vector} or \code{matrix} with variables in columns.}
  \item{block}{\code{vector} or \code{factor} grouping/blocking variable. 
               Must Have length equal to the sample size. 
               All group sample sizes must be greater than 1.}
}
\details{
    If \code{block} describes a unique sample group, ordinary sample standard deviation will be computed.
}
\value{
    \code{numeric} \code{vector} of size the dimensionality where each entry is the variable-wise pooled group sample standard deviation.
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
    Internal function not to be called by end-user.
}
\seealso{
    \itemize{
        \item\code{\link[MVR]{pooled.mean}} Pooled group sample mean
        \item\code{\link[MVR]{mvr}} Mean-Variance Regularization and Variance Stabilization.
        \item\code{\link[MVR]{MeanVarReg}} Mean-Variance Regularization Core Function.
    }
}
\keyword{multiple sample groups}
\keyword{pooled variance estimator}
