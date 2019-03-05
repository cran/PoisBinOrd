\name{gen.PoisBinOrd}
\alias{gen.PoisBinOrd}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Simulates a sample of size n from a set of multivariate Poisson, binary, and
ordinal variables}

\description{This function simulates a sample of size n from a set of multivariate Poisson, binary, and
ordinal data with pre-specified marginals and a correlation matrix.}

\usage{gen.PoisBinOrd(n, n.P, n.B, n.O, lambda.vec = NULL, prop.vec = NULL, 
prop.list = NULL, final.corr.mat)}

\arguments{
  \item{n}{Number of variates.}
  \item{n.P}{Number of Poisson variables.}
  \item{n.B}{Number of binary variables.}
  \item{n.O}{Number of ordinal variables.}
  \item{lambda.vec}{Rate vector for Poisson variables.}
  \item{prop.vec}{Probability vector for binary variables.}
  \item{prop.list}{A list of probability vectors for ordinal variables.}
  \item{final.corr.mat}{Final correlation matrix produced from \code{\link{overall.corr.mat}}}
}


\value{A matrix of size n*(n.P + n.B + n.O), of which the first n.P 
columns are Poisson variables, the next n.B columns are binary variables, and the last n.O 
columns are ordinal variables.}

\examples{
\dontrun{
n=10000
n.P<-2
n.B<-2
n.O<-2
lambda.vec<-sample(10,2)
prop.vec<-runif(2)
prop.list<-list(c(0.3,0.6,0.7),c(0.2,0.3,0.5))
M=c(-0.05, 0.26, 0.14, 0.09, 0.14, 0.12, 0.13, -0.02, 0.17, 0.11,
-0.04, 0.19, 0.10, 0.35, 0.39)
N=diag(6)
N[lower.tri(N)]=M
corr.mat=N+t(N)
diag(corr.mat)<-1
final.corr.mat=overall.corr.mat(n.P,n.B,n.O,lambda.vec,prop.vec,
prop.list,corr.vec=NULL,corr.mat)
mymixdata=gen.PoisBinOrd(n,n.P,n.B,n.O,lambda.vec,prop.vec,prop.list, 
final.corr.mat)
}
}

