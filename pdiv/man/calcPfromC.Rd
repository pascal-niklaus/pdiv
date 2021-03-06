% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trans.R
\name{calcPfromC}
\alias{calcPfromC}
\alias{calcPfromCfun}
\title{Calculate transition matrix from competition matrix}
\usage{
calcPfromC(C)

calcPfromCfun(m, digits = 14)
}
\arguments{
\item{C}{square C matrix}

\item{m}{dimension of C matrix (i.e. number of species)}

\item{digits}{number of digits for coefficients in the resulting code (default 14)}
}
\value{
C matrix, or function yielding C matrix
}
\description{
Determine transition matrix P from competition matrix C using the procedure
described in Ulrich et al. (2014).
}
\details{
Two functions are provided. The first function (\code{calcPfromC})
evaluates the data in C recursively, directly producing matrix
P. Due to the recursion, it is slow in particular for large
matrices.  The second function (\code{calcPfromCcode}) is a function generator that
produces code that explicitly evaluates all cells in matrix P with hard-coded expressions.
Calling this code later is much faster.
}
\examples{
library(pdiv)
## example from Fig. 2 in Ulrich et al.
C4<-matrix(c(1,.8,.3,.2,1,.8,.7,.2,1),byrow=TRUE,nrow=3)
calcPfromC(C4)
#      [,1] [,2] [,3]
# [1,] 0.24 0.72 0.18
# [2,] 0.13 0.16 0.68
# [3,] 0.63 0.12 0.14
fun <- calcPfromCfun(3)
fun(C4)
#      [,1] [,2] [,3]
# [1,] 0.24 0.72 0.18
# [2,] 0.13 0.16 0.68
# [3,] 0.63 0.12 0.14
fun
# function(C)
# {
#   P <- matrix(nrow=nrow(C), ncol=ncol(C))
#   P[1,1] <- 1*C[1,2]*C[1,3]
#   P[1,2] <- 0.5*C[1,2]+0.5*C[1,2]*C[2,3]
#   P[1,3] <- 0.5*C[1,3]+0.5*C[1,3]*C[3,2]
#   P[2,1] <- 0.5*C[2,1]+0.5*C[1,3]*C[2,1]
#   P[2,2] <- 1*C[2,1]*C[2,3]
#   P[2,3] <- 0.5*C[2,3]*C[3,1]+0.5*C[2,3]
#   P[3,1] <- 0.5*C[1,2]*C[3,1]+0.5*C[3,1]
#   P[3,2] <- 0.5*C[2,1]*C[3,2]+0.5*C[3,2]
#   P[3,3] <- 1*C[3,1]*C[3,2]
#   P
# }
}
\references{
Ulrich et al. (2014) Matrix models for quantifying
    competitive intransitivity from species abundance data. Oikos
    123:1057-1070.
}
\author{
Pascal Niklaus \email{pascal.niklaus@ieu.uzh.ch}
}
