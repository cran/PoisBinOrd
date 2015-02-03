\name{intermediate.corr.BO}
\alias{intermediate.corr.BO}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Computes an intermediate normal correlation matrix for any combination of binary and ordinal variables given the specified correlation matrix}

\description{Computes an intermediate normal correlation matrix for any combination of binary and ordinal variables before dichotomization/ordinalization 
given the specified correlation matrix as formulated in Ferrari and Barberio (2012).}

\usage{intermediate.corr.BO(n.B, n.O, prop.vec = NULL, prop.list = NULL, corr.vec = NULL, 
corr.mat = NULL)}

\arguments{
  \item{n.B}{Number of binary variables.}
  \item{n.O}{Number of ordinal variables.}
   \item{prop.vec}{Probability vector for binary variables.}
  \item{prop.list}{A list of probability vectors for ordinal variables.}
  \item{corr.vec}{Vector of elements below the diagonal of correlation matrix ordered column-wise.}
  \item{corr.mat}{Specified correlation matrix.}
}

\value{A correlation matrix of size (n.B+n.O)*(n.B+n.O).}

\references{
Ferrari, P.A. and Barberio, A. (2012). Simulating ordinal data. Multivariate Behavioral Research, 47(4), 566-589.

Barberio, A. and Ferrari, P.A. (2014). GenOrd: Simulation of ordinal and discrete variables with given correlation matrix and marginal distributions.
http://www.cran.r-project.org/web/packages/GenOrd.
}

\seealso{\code{\link{intermediate.corr.P_BO}}}

\examples{
\dontrun{
n.B<-1
n.O<-2
prop.vec=0.7
prop.list = list(cumsum(c(0.30, 0.40)), cumsum(c(0.4, 0.2, 0.3)))
corr.mat = matrix ( c(
1.0000000, 0.1767231, 0.3006186,
0.1767231, 1.0000000, -0.139923,
0.3006186, -0.1399230, 1.0000000),3,3)
intmatBO=intermediate.corr.BO(n.B, n.O, prop.vec, prop.list, corr.vec = NULL, 
corr.mat)

n.B<-1
n.O<-1
prop.vec<-c(0.3)
prop.list<-list(c(0.3,0.6))
corr.mat=matrix(c(1,0.2,0.1,0.2,1,0.5,0.1,0.5,1),3,3)
intmatBO=intermediate.corr.BO(n.B, n.O, prop.vec, prop.list, corr.vec = NULL, 
corr.mat)

n.B<-2
prop.vec=c(0.4,0.7)
corr.mat = matrix(c(1,-0.3,-0.3,1),2,2)
intmatBB=intermediate.corr.BO(n.B, n.O=0, prop.vec, prop.list= NULL, corr.vec = NULL, 
corr.mat)

#See Tetra.Corr.BB in R package BinNonNor 
#Tetra.Corr.BB(n.BB=2, prop.vec=c(0.4,0.7), corr.vec = NULL, corr.mat = corr.mat)

n.B<-0
n.O<-2
prop.list = list(cumsum(c(0.30, 0.40)), cumsum(c(0.4, 0.2, 0.3)))
corr.mat = matrix(c(1.0000000, -0.139923,-0.139923,1.0000000),2,2)
intmatOO=intermediate.corr.BO(n.B, n.O, prop.vec=NULL, prop.list, corr.vec = NULL, 
corr.mat)

#See IntermediateOO(plist, OOCorrMat) in R package OrdNor
#IntermediateOO(plist=list(cumsum(c(0.30, 0.40)), cumsum(c(0.4, 0.2, 0.3))), 
OOCorrMat=corr.mat)
}
}

