\name{is.valid}
\alias{is.valid}
\title{Validation Subroutine to Test the Bootstrap Set of Indices}
\description{
    Internal subroutine called by end-user \code{\link[MVR]{mvrt.test}} function to validate the bootstrap set of indices 
    if a bootstrap test is undertaken in \code{\link[MVR]{mvrt.test}}. Check that the bootstrapped indices of samples 
    (from all-groups samples) are unique per sample group in a multiple sample group situation.
}
\usage{
    is.valid(x, def, ng)
}
\arguments{
  \item{x}{\code{vector} of bootstrap set of indices.}
  \item{def}{\code{list} of samples indices per sample group in original data.}
  \item{ng}{Positive \code{integer} scalar of the number of sample groups.}
}
\details{
    None
}
\value{
    \item{ok}{\code{logical} scalar, taking on value \code{TRUE} if valid and \code{FALSE} otherwise.}
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
    \code{\link[MVR]{mvrt.test}} Mean-Variance Regularized t-Test Statistic with Significance.
}
\keyword{bootstrap}
\keyword{resampling with replacement}


