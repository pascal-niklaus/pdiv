#' Hill's series of diversity indices
#'
#' Hill proposed a unified notation of diversity numbers of different
#' orders \code{a}. As special cases, it also contains the number of
#' species (\code{a=0}), Shannon's index of diversity (\code{a=1}),
#' and the inverse of Simpson's index (\code{a=2}). \code{a=-Inf} and
#' \code{a=+Inf} correspond to the inverse of the proportional
#' abundance of the rarest and the most common species.
#'
#' \code{hill.Na(x,a)} calculated these diversity numbers, whereas
#' \code{hill.Ha(x,a)} computes the generalized entropy measures of
#' order \code{a}.  For convenience, \code{simpson.D} and
#' \code{shannon.H} define these special cases and essentially yield
#' \code{1/hill.Na(x,2)} and \code{hill.Ha(x,1)}.
#'
#' Many measures of diversity indices have been proposed, which have
#' different mathematical properties. A particularly useful one is
#' implemented as \code{evenness.1overD}, and returns \code{1/(DS)},
#' where D is Simpson's index of dominance and S is the total number
#' of species present. It has many desirable properties, including
#' independence from species richness, which is not the case for many
#' other metrics (see Smith and Wilson for a detailed illustration,
#' and example below).
#' 
#' @param x Proportional measure of species abundances. These need not
#' be standardized to unity sum.
#' @param a Order of diversity number to be calculated. Defaults to \code{1}.
#' @return Diversity metric
#' 
#' @references Hill MO (1973) Diversity and evenness: a unifying
#' notation and its consequences. Ecology 54(2):428--432.
#'
#' Smith B & Wilson JB (1996) A consumer's guide to
#' diversity indices, OIKOS 76: 70--82.
#' 
#' @examples
#' library(pdiv)
#' 
#' x1 <- c(1,1,1,3)
#' x2 <- rep(x1,2)
#' 
#' shannon.H(x1)
#' ## [1] 1.242453
#' 
#' simpson.D(x1)
#' ## [1] 0.3333333
#' 
#' hill.Na(x1)
#' ## [1] 3.464102
#' hill.Ha(x1,1)
#' ## [1] 1.242453               # equal to Shannon's H
#' 
#' hill.Na(x1,0)
#' ## [1] 4                      # species numbers
#'
#' hill.Na(x1,2)
#' ## [1] 3                      # equal to 1/D
#'
#' evenness.1overD(x1)           # note that E(1/D) is 
#' ## [1] 0.75                   # equal for x1 and x2
#' evenness.1overD(x2)
#' ## [1] 0.75
#' 
#' shannon.H(x1)/log(length(x1)) # whereas this is 
#' ## [1] 0.8962406              # not the case for
#' shannon.H(x2)/log(length(x2)) # the Shannon-based
#' ## [1] 0.9308271              # evenness index
#' 
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' 
#' @rdname hill
#' @export    
hill.Na <- function(x,a=1)
{
    if(is.vector(x))
        x <- matrix(x,nrow=1)
    sapply(1:nrow(x),
           function(i) {
               .hill.Na(x[i,],a);
           })
}


.hill.Na <- function(x, a=1)
{
    x <- as.numeric(x)
    x <- x[ is.finite(x) & x > 0 ]
    p <- x/sum(x)
    if(a==1) 
        exp(-sum(p*log(p)))
    else 
        sum(p^a)^(1/(1-a))    
}

#' @rdname hill
#' @export
hill.Ha <- function(x, a=1)
{
    log(hill.Na(x,a))
}

#' @rdname hill
#' @export
shannon.H <- function(x)
{
    hill.Ha(x,1);
}

#' @rdname hill
#' @export
simpson.D <- function(x)
{
    1 / hill.Na(x,2)
}

#' @rdname hill
#' @export
evenness.1overD <- function(x)
{
    1 / simpson.D(x) / hill.Na(x,0)
}
