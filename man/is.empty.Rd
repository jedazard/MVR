\name{is.empty}
\alias{is.empty}
\title{Checks if Object is Empty}
\description{
    Internal subroutine called by end-user \code{\link[MVR]{mvr}} and \code{\link[MVR]{mvrt.test}} functions to represent the empty 
    array, matrix, or vector of zero dimension or length. Often returned by expressions and functions whose value is undefined.
}
\usage{
    is.empty(x)
}
\arguments{
  \item{x}{\code{array}, \code{matrix} or \code{vector} of any type.}
}
\details{
    None
}
\value{
    \code{logical} scalar. Returns \code{TRUE} if its argument is empty and \code{FALSE} otherwise.
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
    \code{is.null} from R package \pkg{base}.
}
\keyword{Empty Object}


