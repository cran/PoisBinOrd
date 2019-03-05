\name{correlation.bound.check}
\alias{correlation.bound.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Checks if the pairwise correlation among variables are within the feasible range}

\description{This function checks if there are range violations among correlation of Poisson-Poisson, 
Poisson-binary, Poisson-ordinal, binary-binary, binary-ordinal, and ordinal-ordinal combinations.}

\usage{correlation.bound.check(n.P, n.B, n.O, lambda.vec = NULL, prop.vec = NULL, 
prop.list = NULL, corr.vec = NULL, corr.mat = NULL)}

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

\value{The function returns TRUE if no specification problem is encountered. Otherwise, it returns an error message.}

\references{
Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper
correlation bounds. The American Statistician, 65(2), 104-109.

Demirtas, H., Hedeker, D., and Mermelstein, R.J. (2012). Simulation of massive public health data
by power polynomials. Statistics in Medicine, 31(27), 3337-3346.
}

\seealso{\code{\link{validation.corr}}, \code{\link{correlation.limits}} }


\examples{
\dontrun{
n.P<-1
n.B<-1
n.O<-1
lambda.vec<-c(1)
prop.vec<-c(0.3)
prop.list<-list(c(0.3,0.6))
corr.mat=matrix(c(1,0.2,0.1,0.2,1,0.5,0.1,0.5,1),3,3)
correlation.bound.check(n.P,n.B,n.O,lambda.vec,prop.vec,prop.list,corr.vec=NULL,
corr.mat)

n.P<-2
n.B<-2
n.O<-2
lambda.vec<-c(1,2)
prop.vec<-c(0.3,0.5)
prop.list<-list(c(0.3,0.6),c(0.5,0.6))
corr.mat=matrix(0.8,6,6)
diag(corr.mat)=1
correlation.bound.check(n.P,n.B,n.O,lambda.vec,prop.vec,prop.list,corr.vec=NULL,
corr.mat)
}
}