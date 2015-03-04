\name{withinsumsq}
\alias{withinsumsq}
\title{Within-Cluster Sum of Squares Distances Subroutine}
\description{
    Wrapping subroutine around internal C subroutine for returning Monte-Carlo replicates 
    of within-cluster sum of squares distances for a given number of clusters and
    under a standard Gaussian reference distribution.
}
\usage{
    withinsumsq(n, p, B, k) 
}
\arguments{
  \item{n}{Positive \code{integer} scalar of the number of rows in the data matrix (where points to cluster are by rows (usually samples)).}
  \item{p}{Positive \code{integer} scalar of the number of columns in the data matrix (where variables are by columns).}
  \item{k}{Positive \code{integer} scalar of the fixed number of clusters.}
  \item{B}{Positive \code{integer} scalar of the number of Monte Carlo replicates of the inner loop of the gap statistic function.}
}
\details{
    None
}
\value{
    \item{lWk.mc}{\code{numeric} \code{B}-vector of Monte Carlo replicates of within-cluster sum of squares.}
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
    Internal subroutine of \code{\link[MVR]{sim.dis}}, not to be called by end-user.
}
\seealso{
    \code{\link[MVR]{sim.dis}} Similarity Statistic.
}
\keyword{cluster}
\keyword{distance}
\keyword{dispersion}
