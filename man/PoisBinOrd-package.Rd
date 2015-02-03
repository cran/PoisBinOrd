\name{PoisBinOrd-package}
\alias{PoisBinOrd-package}
\alias{PoisBinOrd}
\docType{package}
\title{Data generation with count, binary and ordinal components}
\description{Provides R functions for generation of multiple count, binary and ordinal variables simultaneously 
             given the marginal characteristics and association structure.}
\details{
\tabular{ll}{
Package: \tab PoisBinOrd\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-03\cr
License: \tab GPL (>= 2) \cr
}PoisBinOrd package consists of ten functions. The functions \code{\link{validation.bin}}, \code{\link{validation.ord}}, and
\code{\link{validation.corr}} validate the specified quantities to avoid users from obvious specification errors. 
\code{\link{correlation.limits}} returns the lower and upper bounds of the pairwise correlation of Poisson-Poisson, Poisson-binary, Poisson-ordinal, binary-binary, binary-ordinal, and ordinal-ordinal combinations given their marginal distributions,
i.e. returns the range of feasible pairwise correlations. The function \code{\link{correlation.bound.check}}
checks the validity of the values of pairwise correlations. The functions \code{\link{intermediate.corr.PP}}, 
\code{\link{intermediate.corr.BO}}, and  \code{\link{intermediate.corr.P_BO}} 
computes intermediate correlation matrix for Poisson-Poisson combinations, binary/ordinal and binary/ordinal combinations, 
and Poisson and binary/ordinal combinations, respectively. The function \code{\link{overall.corr.mat}} assembles 
the final correlation matrix. The engine function \code{\link{gen.PoisBinOrd}}
generates mixed data in accordance with the specified marginal and correlational quantities.
}

\author{Gul Inan, Hakan Demirtas 

Maintainer: Gul Inan <ginan@metu.edu.tr>}

\references{Amatya, A. and Demirtas, H. (2015). Simultaneous generation of multivariate mixed data with Poisson and normal
marginals. Journal of Statistical Computation and Simulation. In press.

Demirtas, H. and Hedeker, D. (2011). A practical way for computing approximate lower and upper
correlation bounds. The American Statistician, 65(2), 104-109.

Demirtas, H., Hedeker, D., and Mermelstein, R.J. (2012). Simulation of massive public health data
by power polynomials. Statistics in Medicine, 31(27), 3337-3346.

Ferrari, P.A. and Barberio, A. (2012). Simulating ordinal data. Multivariate Behavioral Research, 47(4), 566-589.
}

\keyword{concurrent generation of Poisson, binary and ordinal variables}
\keyword{generating multivariate Poisson variables}
\keyword{generating multivariate binary variables}
\keyword{generating multivariate ordinal variables}
