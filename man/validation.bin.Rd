\name{validation.bin}
\alias{validation.bin}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Validates the marginal specification of the binary variables}

\description{Checks whether the marginal specification of the binary part is valid and consistent.}

\usage{validation.bin(n.B, prop.vec = NULL)}

\arguments{
  \item{n.B}{Number of binary variables.}
  \item{prop.vec}{Probability vector for binary variables.}
}

\value{The function returns TRUE if no specification problem is encountered. Otherwise, it returns an error message.}

\examples{
n.B<-3
prop.vec<-c(0.25,0.5,0.75)
validation.bin(n.B, prop.vec)

\dontrun{
n.B<-3
validation.bin(n.B)

n.B<-NULL
prop.vec<-c(0.25,0.5,0.75)
validation.bin(prop.vec)

n.B<--3
prop.vec<-c(0.25,0.5,0.75)
validation.bin(n.B, prop.vec)

n.B<-0
prop.vec<-c(0.25,0.5,0.75)
validation.bin(n.B, prop.vec)

n.B<-5
prop.vec<-c(0.25,0.5,0.75)
validation.bin(n.B, prop.vec)

n.B<-3
prop.vec<-c(0.25,0.5,-0.75)
validation.bin(n.B, prop.vec)
}
}
