#' Rarefy species diversity
#'
#' Given a vector of abundances, \code{rarefy} and \code{ram_rarefy}
#' calculate species richness rarefied to a smaller sample size.
#' \code{rarefy.sim} is an alternative function that uses repeated
#' random subsampling. It allows to compute a used-defined function
#' for each subsample, instead of just species richness.
#'
#' \code{rarefy} is fast and uses a logarithmic representation of
#' factorials to avoid overflow issues even for very large samples as
#' they are typical in e.g. the analysis of large sequence data.
#'
#' \code{rarefy2} does the same but uses function \code{lchoose2}
#' instead of the standard \code{lchoose} to calculate binomial
#' coefficient. \code{rarefy2} has the advantage that it works for
#' fractional sample sizes, which can be useful for
#' interpolation. \code{lchoose2} is similar to \code{lchose} but
#' accepts fractional numbers. \code{rarefy2} is slightly slower than
#' \code{rarefy}.
#'
#' \code{rarefy.sim} takes random subsamples of the requested size. By
#' default, it returns species richness and thus provides no advantage
#' over \code{rarefy}. However, a user defined function \code{fun} can
#' be passed, which is applied to the random subsample. This allows,
#' for example, to compute the Shannon diversity index of a rarefied
#' sample.
#'
#' \code{rarefy.sim} does not adopt parallel processing
#' because this is more efficiently done at higher levels
#' (e.g. parallelize the computation of rarefied richness of multiple
#' samples).
#'
#' @param abu vector of abundances of the different species (or
#'     OTUs). These must be positive integers; zeroes are discarded.
#'
#' @param size size of downsampled population. Note that \code{size}
#'     must be smaller or equal to \code{sum(abu)}, otherwise a
#'     warning is issued and the result is NA. For \code{rarefy} (but
#'     not \code{rarefy.sim}), size may be a vector.
#'
#' @param n,k arguments to determine bionomial coefficient, can be fractional
#'
#' @param se logical indicating whether a standard error of rarefied
#'     richness should be returned (only \code{rarefy.sim}). Note that
#'     this is only the standard error of the simulation process, and
#'     does not factor in that the original abundances also are a
#'     sample subject to uncertainty!
#'
#' @param fun function to be applied to random subsample. As argument,
#'     it must take a vector with abundances of the species remaining
#'     in the subsample. The optional arguments \code{...} are also
#'     passed.
#'
#' @param ... optional arguments passed down to \code{fun}.
#'
#' @param nsim number of random samples that are drawn (only
#'     \code{rarefy.sim}).
#'
#' @return rarefied richness, or a list containing rarefied richness
#'     and its standard error (only \code{rarefy.sim}, and only if
#'     \code{se=TRUE}).
#'
#' @examples
#' library(pdiv)
#' abu <- c(10,100,1000,10000)
#' rarefy(abu,10)
#' rarefy.sim(abu, 10, nsim=1000)
#'
#' sizes <- seq(0, 100, by=2)
#' rich <- rarefy(abu, sizes)
#' plot(rich ~ sizes)
#'
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' @importFrom stats sd
#' @rdname rarefy
#' @export
rarefy <- function(abu, size) {
    if(!all(is.finite(abu)))
        return(NA)
    abu <- abu[abu > 0]  # remove zero abundance species
    K <- length(abu)     # number of species/OTUs
    N <- sum(abu)        # number of individuals
    sapply(size,
           function(size) {
               if(is.na(size) || size < 0  || size > N) {
                   warning("Rarefaction must reduce sample size, returning NA!")
                   return(NA)
               }
               r <- sapply(1:K, function(i) lchoose(N-abu[i], size))
               K-sum(exp(r-lchoose(N,size)))
           })
}

#' @rdname rarefy
#' @export
rarefy2 <- function(abu, size)
{
    if(!all(is.finite(abu)))
        return(NA)
    abu <- abu[abu > 0]
    K <- length(abu)
    N <- sum(abu)
    sapply(size,
           function(size) {
               if(is.na(size) || size < 0  || size > N) {
                   warning("Rarefaction must reduce sample size, returning NA!")
                   return(NA)
               }
               logk <- lchoose2(N, size)
               K- sum(sapply(1:K,
                             function(i) exp(lchoose2(N-abu[i],size)-logk)))
           })
}

#' @rdname rarefy
#' @export
lchoose2 <- function(n, k) {
    sum(lfactorial(c(n, k, n-k)) * c(1, -1, -1))
}

#' @rdname rarefy
#' @export
rarefy.sim <- function(abu, size, nsim=1000, se=TRUE, fun=NULL, ...) {
    if(!all(is.finite(abu)))
        return(NA)
    abu <- abu[abu > 0]  # remove zero abundance species
    N <- sum(abu)       # individuals
    if(is.na(size) || size < 0  || size > N) {
        warning("Rarefaction must reduce sample size, returning NA!")
        return(NA)
    }
    if(size == N)
        return(list(S=length(abu), se=0))
    breaks <- cumsum(sort(abu, decreasing=TRUE))+1
    ivec <- as.integer(1:N)
    r <- sapply(1:nsim,
                function(i) {
                    s <- sample(ivec, size)
                    tmp <- findInterval(s,breaks)
                    if(is.null(fun))
                        length(unique(tmp))
                    else
                        fun(sort(as.numeric(table(tmp))), ...)
                })
    if(se)
        list(S=mean(r), se=sd(r)/sqrt(nsim))
    else
        mean(r)
}
