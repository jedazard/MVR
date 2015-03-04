\name{merging.cluster}
\alias{merging.cluster}
\title{Variable-Cluster Configuration Merging Subroutine}
\description{
    Merging subroutine of variable-cluster configurations, determined for each group.
    Internal function called by \code{\link[MVR]{mvr}} and internally by \code{\link[MVR]{mvrt.test}}. 
    Takes as argument the cluster membership matrix for each variable by group, 
    which is returned by the membership value of \code{\link[MVR]{MeanVarReg}}.
}
\usage{
    merging.cluster(M)
}
\arguments{
  \item{M}{Cluster membership matrix for each variable by group, 
           where variables are by rows, and groups by columns.
           Returned by the membership value of \code{\link[MVR]{MeanVarReg}}.}
}
\details{
    None
}
\value{
    \item{clus}{\code{numeric} \code{vector} representing the merged-clusters membership of each variable.}
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
        \item\code{\link[MVR]{mvr}} Mean-Variance Regularization and Variance Stabilization.
        \item\code{\link[MVR]{mvrt.test}} Mean-Variance Regularized t-Test Statistic with Significance.
        \item\code{kmeans} (R package \pkg{stats}) for partition clustering in R.
    }
}
\keyword{clustering}
\keyword{membership}

