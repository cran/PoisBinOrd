\name{overall.corr.mat}
\alias{overall.corr.mat}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Computes the final intermediate correlation matrix}

\description{This function computes the final correlation matrix by combining pairwise intermediate correlation matrix entries for Poisson-Poisson, 
Poisson-binary, Poisson-ordinal, binary-binary, binary-ordinal, and ordinal-ordinal combinations. If
the resulting correlation matrix is not positive definite, a nearest positive matrix will be used.}

\usage{overall.corr.mat(n.P, n.B, n.O, lambda.vec = NULL, prop.vec = NULL, prop.list = NULL, 
corr.vec = NULL, corr.mat = NULL)}

\arguments{
  \item{n.P}{Number of Poisson variables.}
  \item{n.B}{Number of binary variables.}
  \item{n.O}{Number of ordinal variables.}
  \item{lambda.vec}{Rate vector for Poisson variables.}
  \item{prop.vec}{Probability vector for binary variables.}
  \item{prop.list}{A list of probability vectors for ordinal variables.}
  \item{corr.vec}{Vector of elements below the diagonal of correlation matrix ordered column-wise.}
  \item{corr.mat}{Specified correlation matrix.}
}


\value{A matrix of size (n.P+n.B+n.O)*(n.P+n.B+n.O)}

\seealso{\code{\link{intermediate.corr.PP}}, \code{\link{intermediate.corr.BO}},\code{\link{intermediate.corr.PBO}}}

\examples{
\dontrun{
n.P<-1
n.B<-1
n.O<-1
lambda.vec<-c(1)
prop.vec<-c(0.3)
prop.list<-list(c(0.3,0.6))
corr.vec=NULL
corr.mat=matrix(c(1,0.2,0.1,0.2,1,0.5,0.1,0.5,1),3,3)
finalmat=overall.corr.mat(n.P,n.B,n.O,lambda.vec,prop.vec,prop.list,corr.vec=NULL, 
corr.mat)
}
}
