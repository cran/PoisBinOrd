\name{intermediate.corr.P_BO}
\alias{intermediate.corr.P_BO}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Computes the pairwise entries of the intermediate normal correlation matrix for all Poisson-binary and Poisson-ordinal variable 
combinations given the specified correlation matrix.}

\description{This function computes the pairwise entries of the intermediate normal correlation matrix for all Poisson-binary and Poisson-ordinal variable 
combinations given the specified correlation matrix as formulated in Amatya and Demirtas (2015). }

\usage{intermediate.corr.P_BO(n.P, n.B, n.O, lambda.vec = NULL, prop.vec = NULL, 
prop.list = NULL, corr.vec = NULL, corr.mat = NULL)}

\arguments{
  \item{n.P}{Number of Poisson variables.}
  \item{n.B}{Number of binary variables.}
  \item{n.O}{Number of ordinal variables.}
  \item{lambda.vec}{Rate vector for Poisson variables}
  \item{prop.vec}{Probability vector for binary variables.}
  \item{prop.list}{A list of probability vectors for ordinal variables.}
  \item{corr.vec}{Vector of elements below the diagonal of correlation matrix ordered column-wise.}
  \item{corr.mat}{Specified correlation matrix.}
}

\value{A matrix of n.P*(n.B+n.O)}

\references{
Amatya, A. and Demirtas, H. (2015). Simultaneous generation of multivariate mixed data with Poisson and normal
marginals. Journal of Statistical Computation and Simulation. In press.

}

\seealso{\code{\link{intermediate.corr.PP}}, \code{\link{intermediate.corr.BO}}}

\examples{
\dontrun{
n.P<-1
n.B<-1
n.O<-1
lambda.vec<-c(1)
prop.vec<-c(0.3)
prop.list<-list(c(0.3,0.6))
corr.mat=matrix(c(1,0.2,0.1,0.2,1,0.5,0.1,0.5,1),3,3)

intmatP_BO=intermediate.corr.P_BO(n.P, n.B, n.O, lambda.vec, prop.vec, prop.list, 
corr.vec = NULL, corr.mat)
}
}