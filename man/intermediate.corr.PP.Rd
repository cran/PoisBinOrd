\name{intermediate.corr.PP}
\alias{intermediate.corr.PP}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Computes an intermediate normal correlation matrix for Poisson variables given the specified correlation matrix}

\description{This function computes the intermediate normal correlation matrix for Poisson-Poisson combinations before inverse cdf matching 
as formulated in Amatya and Demirtas (2015).}

\usage{intermediate.corr.PP(n.P, lambda.vec, corr.vec = NULL, corr.mat = NULL)}

\arguments{
 \item{n.P}{Number of Poisson variables}
 \item{lambda.vec}{Rate vector for Poisson variables}
 \item{corr.vec}{Vector of elements below the diagonal of correlation matrix ordered column-wise.}
 \item{corr.mat}{Specified correlation matrix.}
}

\value{A correlation matrix of size n.P*n.P}

\references{Amatya, A. and Demirtas, H. (2015). Simultaneous generation of multivariate mixed data with Poisson and normal
marginals. Journal of Statistical Computation and Simulation, 85(15), 3129-3139.
}

\seealso{\code{\link{intermediate.corr.P_BO}}}


\examples{
n.P<-3
lambda.vec<-c(1,2,3)
corr.mat<-matrix(c(1,0.352,0.265,0.352,1,0.121,0.265,0.121,1),n.P,n.P)
intmatPP=intermediate.corr.PP(n.P,lambda.vec,corr.vec=NULL,corr.mat)

\dontrun{
#See also cmat.star in  R package PoisNor 
#cmat.star(no.pois=3,no.norm=0,corMat=corr.mat,lamvec=lambda.vec)
}
}
