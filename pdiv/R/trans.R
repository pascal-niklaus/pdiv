#' Calculate transitivity from transition matrix
#'
#' Given a transition matrix P, this function calculates transitivity tau_P as defined
#' in equation 19 in Ulrich et al. (2014).
#'
#' @param P square P matrix
#' @return transitivity measure `tau_P'
#' 
#' @references Ulrich et al. (2014) Matrix models for quantifying
#'     competitive intransitivity from species abundance data. Oikos
#'     123:1057-1070.
#'
#' @examples
#' library(pdiv)
#' ## example from Fig. 2 in Ulrich et al.
#' C1<-matrix(c(1,1,1,0,1,1,0,0,1),byrow=TRUE,nrow=3)
#' C4<-matrix(c(1,.8,.3,.2,1,.8,.7,.2,1),byrow=TRUE,nrow=3)
#' P1<-calcPfromC(C1)
#' P4<-calcPfromC(C4)
#' calcTp(P1)
#' calcTp(P4) 
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' @export  
calcTp <- function(P)
{
    m <- nrow(P)
    n <- 0
    for(i in 1:m) 
        for(j in 1:m) 
            for(k in (i+1):m) 
                if(!(i>=k || i==j || k==j || k>m)) 
                    if(P[i,j] <- P[k,j]) n <- n + 1                
    1-2*n/(m*(m-1)*(m-2))
}

#' Calculate transition matrix from competition matrix
#'
#' Determine transition matrix P from competition matrix C using the procedure
#' described in Ulrich et al. (2014).
#'
#' Two functions are provided. The first function (\code{calcPfromC})
#' evaluates the data in C recursively, directly producing matrix
#' P. Due to the recursion, it is slow in particular for large
#' matrices.  The second function (\code{calcPfromCcode}) is a function generator that
#' produces code that explicitly evaluates all cells in matrix P with hard-coded expressions.
#' Calling this code later is much faster.
#' 
#' @param C square C matrix
#' @param m dimension of C matrix (i.e. number of species)
#' @param digits number of digits for coefficients in the resulting code (default 14)
#' @return C matrix, or function yielding C matrix
#' 
#' @references Ulrich et al. (2014) Matrix models for quantifying
#'     competitive intransitivity from species abundance data. Oikos
#'     123:1057-1070.
#'
#' @examples
#' library(pdiv)
#' ## example from Fig. 2 in Ulrich et al.
#' C4<-matrix(c(1,.8,.3,.2,1,.8,.7,.2,1),byrow=TRUE,nrow=3)
#' calcPfromC(C4)
#' #      [,1] [,2] [,3]
#' # [1,] 0.24 0.72 0.18
#' # [2,] 0.13 0.16 0.68
#' # [3,] 0.63 0.12 0.14
#' fun <- calcPfromCfun(3)
#' fun(C4)
#' #      [,1] [,2] [,3]
#' # [1,] 0.24 0.72 0.18
#' # [2,] 0.13 0.16 0.68
#' # [3,] 0.63 0.12 0.14
#' fun
#' # function(C)
#' # {
#' #   P <- matrix(nrow=nrow(C), ncol=ncol(C))
#' #   P[1,1] <- 1*C[1,2]*C[1,3]
#' #   P[1,2] <- 0.5*C[1,2]+0.5*C[1,2]*C[2,3]
#' #   P[1,3] <- 0.5*C[1,3]+0.5*C[1,3]*C[3,2]
#' #   P[2,1] <- 0.5*C[2,1]+0.5*C[1,3]*C[2,1]
#' #   P[2,2] <- 1*C[2,1]*C[2,3]
#' #   P[2,3] <- 0.5*C[2,3]*C[3,1]+0.5*C[2,3]
#' #   P[3,1] <- 0.5*C[1,2]*C[3,1]+0.5*C[3,1]
#' #   P[3,2] <- 0.5*C[2,1]*C[3,2]+0.5*C[3,2]
#' #   P[3,3] <- 1*C[3,1]*C[3,2]
#' #   P
#' # }
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#'
#' @rdname calcPfromC
#' @export  
calcPfromC <- function(C)
{
    m <- nrow(C)
    P <- matrix(nrow=m, ncol=m)
    
    for(i in 1:m) 
        for(j in 1:m)
            P[i,j] <- walkPermutations(i,j,C)        
    P
}

walkPermutations <- function(end, orig, C)
{
    walkPermutationsR(setdiff(1:ncol(C),orig),end,orig,C)
}

walkPermutationsR <- function(rset, end, orig, C)
{
    n <- length(rset)
    x <- sapply(rset,
                function(i) {
                    if(i == end) {
                        1/n * C[end,orig]
                    } else if(n == 1) {
                        1/n * C[orig,i]
                    } else {
                        1/n * C[orig,i] * walkPermutationsR(
                                              setdiff(rset,i),
                                              end,orig,
                                              C)
                    }
                })
    sum(x)
}

#' @rdname calcPfromC
#' @export  
calcPfromCfun <- function(m,digits=14)
{
    body <- paste(
        sapply(1:m, function(i)
            paste(
                sapply(1:m,
                       function(j)
                           sprintf("  P[%d,%d] <- %s",i,j,walkPermutationsCode(i,j,m,digits))),
                collapse="\n")),
        collapse="\n")
    txt <- sprintf("function(C)\n{\n  P <- matrix(nrow=nrow(C), ncol=ncol(C))\n%s\n  P\n}\n\n",body)
    return(eval(parse(text=txt)))
}

walkPermutationsCode <- function(end,orig,m,digits=14)
{
    r <- walker(setdiff(1:m,orig), end, orig)
    paste(
        sapply(seq_along(r),
               function(i)
                   sprintf("%.*g*%s",digits,r[[i]],names(r)[[i]])),
        collapse="+")
}

addTerm <- function(coefList, terms, coef) {
    key <- paste(sort(terms),collapse="*")
    coefList[[key]] <- c(coefList[[key]],0)[1] + coef
    coefList
}

walker <- function(rset,end,orig,terms=c(),f=1.0,coefList=list()) {
    n <- length(rset)
    for(i in rset) {
        if(i == end)
            coefList <- addTerm(coefList,c(terms,sprintf("C[%d,%d]",end,orig)),f/n) 
        else if(n == 1)
            coefList <- addTerm(coefList,c(terms,sprintf("C[%d,%d]",orig,i)),f/n) 
        else
            coefList <- walker(rset=setdiff(rset,i),
                               end=end,
                               orig=orig,
                               terms=c(terms,sprintf("C[%d,%d]",orig,i)),
                               f=f/n,
                               coefList=coefList)
    }
    coefList        
}
