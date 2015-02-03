\name{validation.ord}
\alias{validation.ord}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Validates the marginal specification of the ordinal variables}

\description{Checks whether the marginal specification of the ordinal part is valid and consistent.}

\usage{validation.ord(n.O, prop.list = NULL)}

\arguments{
  \item{n.O}{Number of ordinal variables.}
  \item{prop.list}{A list of probability vectors corresponding to each ordinal variable. The i-th element of prop.list is a vector of the cumulative probabilities defining the marginal
distribution of the i-th ordinal component of the multivariate variables. If the i-th ordinal variable has k categories, the i-th vector of the prop.list will contain
k-1 probability values. The k-th element is implicitly 1.}
}

\value{The function returns TRUE if no specification problem is encountered. Otherwise, it returns an error message.}

\examples{
n.O<-3
prop.list<-list(c(0.3,0.6),c(0.25,0.5,0.75),c(0.1,0.2,0.8,0.9))
validation.ord(n.O, prop.list)

\dontrun{
n.O<-3
validation.ord(n.O)

n.O<-NULL
prop.list<-list(c(0.3,0.6),c(0.25,0.5,0.75),c(0.1,0.2,0.8,0.9))
validation.ord(prop.list=prop.list)

n.O<--3
prop.list<-list(c(0.3,0.6),c(0.25,0.5,0.75),c(0.1,0.2,0.8,0.9))
validation.ord(-3, prop.list)

n.O<--0
prop.list<-list(c(0.3,0.6),c(0.25,0.5,0.75),c(0.1,0.2,0.8,0.9))
validation.ord(n.O, prop.list)

n.O<-5
prop.list<-list(c(0.3,0.6),c(0.25,0.5,0.75),c(0.1,0.2,0.8,0.9))
validation.ord(n.O, prop.list)

n.O<-3
prop.list<-list(c(0.3,0.6),c(0.25,0.5,-0.75),c(0.1,0.2,0.8,1.5))
validation.ord(n.O, prop.list)

n.O<-3
prop.list<-list(0.3,c(0.3,0.4),c(0.4,0.2,0.3))
validation.ord(n.O, prop.list)
}
}

