#' Calculate phylogenetic diversity as mean phylogenetic distance of pairs
#'
#' Given a phylogenetic tree and data on community composition, the
#' average distance between all species pairs within a community is
#' calculated.  This procedure optionally includes abundance-weighing
#' (weights are the product of the abundances of the species of each pair).
#' The distance metric used is the sum of the branch lengths from the
#' two species' tips up to the last common ancestor, multiplied by a
#' user-defined scale factor.
#'
#' On platforms other than Windows, \code{calcMPD} parallelizes the
#' time-critical steps if the \code{parallel} package is
#' available. Parallelization is not available on Windows because of a
#' lack of \code{fork} functionality on which \code{mclapply} bases.
#' 
#' @param communities data frame containing a column specifying the
#' communities (e.g. plots or other groups to compare), a column with
#' species names, and optionally a column with species abundances.
#' Species names must match the tip names in the phylogenetic tree.
#' The abundances column is only used if \code{abundance = TRUE}.
#' @param format indicates in which format community composition is
#' passed. Can be either \code{default} or \code{matrix}.
#' @param tree phylogenetic tree of class \code{phylo}
#' @param which indicates which results to return
#' (\code{default="MPD"}). The vector may contain 
#' \code{"MPD"}, and \code{"dist"}. \code{"MPD"}
#' indicates mean phylogenetic diversity (i.e. mean age of last common ancestor),
#' \code{"dist"} will return the age of the last common ancestor for all pairs.
#' @param scale a user-defined scale factor (default 0.5) by which the
#' total branch length of the connection between a pair is
#' scaled. This is useful in ultrametric trees to arrive at the
#' distance to the last common ancestor since the branch length from
#' both tips would be cumulated otherwise.
#' @param abundance logical indicating whether abundances are used
#' to weigh the distance between each pair.
#' @return a list containing the items selected by \code{which}.
#' @seealso \code{\link{calcPD}}
#' @seealso \code{\link{communityToMatrix}}
#' @seealso \code{\link{communityFromMatrix}}
#' @references Warwick RM and Clark KR (1995) New 'biodiversity' measures
#' reveal a decrease in taxonomic distinctness with increasing stress.
#' Mar Ecol Prog Ser 29:301-305.
#' @examples
#'
#' library(pdiv)
#' data(composition.example)
#' data(phylotree.example)
#'
#' plot(phylotree.example)
#'
#' calcMPD(composition.example,phylotree.example)
#' ## $MPD
#' ## [1] 333  85 293 308   0   0   0   0
#'
#' ex2 <- composition.example
#' ex2$abu<-runif(nrow(composition.example))
#' calcMPD(ex2,phylotree.example,which="MPD",abundance=TRUE)
#' 
#' calcMPD(composition.example,phylotree.example,which="dist")
#' ## $dist
#' ## a|j a|b a|c a|d a|e a|f a|h b|c b|d b|e b|f b|h c|d
#' ## 385 295 295 385 385 385 385 250 385 385 385 385 385
#' ## f|j f|h h|j c|e c|f c|h d|e d|h e|j e|h
#' ## 109 109  37 385 385 385  97 109 109 109
#' 
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' @importFrom utils combn
#' @importFrom stats weighted.mean
#' @export 
calcMPD <- function(
    communities = NULL,
    tree = NULL,
    which = c("MPD"),
    format = "default",
    scale = 0.5,
    abundance = FALSE) 
{
    ## Windows can't fork and does not support mclapply
    parallel <- ( Sys.info()["sysname"] != "Windows" ) && requireNamespace("parallel")
    n.cores <- if(parallel) parallel::detectCores() else 1;
    mypar <-
        if(parallel) {
            function(X, FUN) { parallel::mclapply(X,FUN,mc.cores=n.cores) }
        } else {
            lapply
        }

    rr<-list();

    which.options <- c("MPD","dist");
    which <- which.options[ pmatch( which, which.options ) ]

    formats <- c("default", "matrix");
    format <- pmatch(format, formats)
    if (is.na(format)) 
        stop("invalid format for community composition")
    if (format == -1) 
        stop("ambiguous format for community composition")
    if (abundance && (ncol(communities)!=3 || any(!is.numeric(communities[,3]))))
        stop("abundances must be given as third column of argument communities")

    communities <- as.data.frame(communities);  
    if(format %in% "matrix") 
        communities <-
            if(abundance) {
                communityFromMatrix(communities, abundance.column.name="abu")
            } else {
                communityFromMatrix(communities)
            }

    all.sp <- sort(unique(as.character(communities[,2])));
    if(!(all(all.sp %in% tree$tip.label))) 
        stop("Species",paste(all.sp[!(all.sp %in% tree$tip.label)],collapse=", "),"missing in tree")

    ## walk up tree from all tips, collect edge indices
    ntip <- length(tree$tip.label);
    
    branches <- if(parallel) {
        parallel::mclapply(1:ntip, function(x) .walk_to_root(tree,x))
    } else {
        lapply(1:ntip,function(x) .walk_to_root(tree,x));
    }
    names(branches) <- tree$tip.label

    r <- data.frame(com=sort(unique(communities[,1])),MPD=NA)

    num_sp_vec <- match(communities[,2],tree$tip.label)    

    ## determine all pairs per plot
    pairs_by_plot<- mypar(
        r$com,
        function(x) {
            sp.set <- sort(unique(num_sp_vec[communities[,1]==x]))  
            if(length(sp.set)<2) {
                c()
            } else {
                apply(combn(sp.set,2),2,function(x) paste(x,collapse=":"))
            }
        })

    ## create parallel structure with weights
    abu_prod <- 
        mypar(
            r$com,
            function(com) {
                idx <- ( communities[,1] == com ) 
                sapply(pairs_by_plot[[com]],
                       function(pair) {
                           if(abundance) {
                               sp <- unlist(strsplit(pair, ':'))
                               prod( communities[, 3][ idx ][ match(sp, num_sp_vec[idx]) ] )
                           } else {
                               1
                           }
                       })
            }
            )

    ## all unique pairs for which we need a distance
    used_pairs <- sort(unique(unlist(pairs_by_plot)))

    ## phylogenetic distance of each pair
    pairdists <- unlist(mypar(
        used_pairs,
        function(pair) {
            x <- as.numeric(as.character(unlist(strsplit(pair,":"))))

            br <- c(branches[[x[1]]],branches[[x[2]]])
            brc <- table(br)
            br_unique <- as.numeric(as.character(names(brc)[brc==1]))
            sum(tree$edge.length[br_unique])*scale              
        }
        ))

    ## calculate mean phylogenetic distance by plot
    if("MPD" %in% which) {
        rr$MPD <- data.frame(com=r$com,MPD=NA)
        rr$MPD$MPD <-
            sapply(
                seq_along(pairs_by_plot),
                function(i) {
                    pairs <- unlist( pairs_by_plot[i] )
                    if(is.null(pairs)) {
                        0
                    } else {                        
                        weighted.mean(
                            unlist(pairdists[used_pairs %in% pairs]),
                            unlist(abu_prod[i])
                            )
                        
                    }                                 
                }
                )
    }
    
    if("dist" %in% which) {
        ## add names to pair distances
        ## -  this is not so nice since it duplicates the strsplit from above 
        names(pairdists) <- sapply(
            used_pairs,
            function(x) paste(
                tree$tip.label[ as.numeric(as.character(unlist(strsplit(x,":")))) ],
                collapse="|")
            )                       

        rr$dist <- pairdists
    }

    rr;
}

