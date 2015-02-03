\name{correlation.limits}
\alias{correlation.limits}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Computes lower and upper correlation bounds for each pair of variables}

\description{This function computes lower and upper limits for pairwise correlations of Poisson-Poisson, 
Poisson-binary, Poisson-ordinal, binary-binary, binary-ordinal, and ordinal-ordinal combinations.}

\usage{correlation.limits(n.P, n.B, n.O, lambda.vec = NULL, prop.vec = NULL, 
prop.list = NULL)}

\arguments{
  \item{n.P}{Number of Poisson variables.}
  \item{n.B}{Number of binary variables.}
  \item{n.O}{Number of ordinal variables.}
  \item{lambda.vec}{Rate vector for Poisson variables}
  \item{prop.vec}{Probability vector for binary variables.}
  \item{prop.list}{A list of probability vectors for ordinal variables.}
}

\details{While the function computes the exact lower and upper bounds for pairwise correlations among
binary-binary variables as formulated in Demirtas et al. (2012), it computes approximate lower and
upper bounds for pairwise correlations among Poisson-Poisson, Poisson-binary, Poisson-ordinal, binary-ordinal, and ordinal-ordinal variables 
through the method suggested by Demirtas and Hedeker (2011).}

\value{The function returns a matrix of size (n.P + n.B + n.O)*(n.P + n.B + n.O), where the lower triangular
part of the matrix contains the lower bounds and the upper triangular part of the matrix contains the
upper bounds of the feasible correlations.}

\references{
Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper
correlation bounds. The American Statistician, 65(2), 104-109.

Demirtas, H., Hedeker, D., and Mermelstein, R.J. (2012). Simulation of massive public health data
by power polynomials. Statistics in Medicine, 31(27), 3337-3346.
}

\seealso{\code{\link{validation.corr}}, \code{\link{correlation.bound.check}} }

\examples{
\dontrun{
n.P<-3
n.B<-2
n.O<-3
lambda.vec<-c(1,2,3)
prop.vec<-c(0.3,0.5)
prop.list<-list(c(0.3,0.6),c(0.25,0.5,0.75),c(0.1,0.2,0.8,0.9))

#Correlation limits among Poisson variables
correlation.limits(n.P,n.B=0,n.O=0, lambda.vec,prop.vec=NULL, prop.list=NULL)

#See also Cor.PP.Limit in R package PoisNor 

#Correlation limits among binary variables
correlation.limits(n.P=0,n.B,n.O=0,lambda.vec=NULL, prop.vec, prop.list=NULL)

#See also correlation.limits in R package BinNonNor 

#Correlation limits among ordinal variables
correlation.limits(n.P=0,n.B=0,n.O,lambda.vec=NULL, prop.vec=NULL, prop.list)

#See also Limit_forOO in R package OrdNor 

#Correlation limits among Poisson and binary variables and within themselves.
correlation.limits(n.P,n.B,n.O=0,lambda.vec, prop.vec, prop.list=NULL)

#Correlation limits among Poisson and ordinal variables and within themselves.
correlation.limits(n.P,n.B=0,n.O,lambda.vec, prop.vec=NULL, prop.list)

#Correlation limits among binary and ordinal variables and within themselves.
correlation.limits(n.P=0,n.B,n.O,lambda.vec=NULL, prop.vec, prop.list)

#Correlation limits among Poisson, binary, and ordinal variables and within themselves.
correlation.limits(n.P,n.B,n.O,lambda.vec, prop.vec, prop.list)

n.P<-2
lambda.vec=c(-1,1)
correlation.limits(n.P,n.B=0,n.O=0, lambda.vec,prop.vec=NULL, prop.list=NULL)
}
}

