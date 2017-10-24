#' Calculate Petchey \& Gaston's functional diversity index
#'
#' Calculates functional diversity sensu Petchey & Gaston, including
#' the changes suggested in their 2006 update.
#'
#' Note that this function converts the tree obtained from the
#' hierarchical cluster analysis (\code{hcluse}) to a tree of class
#' \code{phylo}, using \code{as.phylo}.  This introduces a scaling
#' factor of two. As a result, the absolute measure \code{FD} is only
#' half the one obtained by Owen Petchey and Jens Schumacher's
#' code. Standardized FD is equal, though.
#' 
#' @param communities data frame containing a column with species
#'     names and a column designating the communities (e.g. plots or
#'     other groups to compare)
#'
#' @param traits data frame containing the traits of each species
#'
#' @param distance the method used to calculate distances in trait
#'     space. Options are \code{"euclidean"}, \code{"maximum"},
#'     \code{"manhattan"}, \code{"canberra"}, \code{"binary"},
#'     \code{"minkowski"}, and \code{"gower"}. See \code{\link{dist}}
#'     and \code{\link{daisy}} for details.
#'
#' @param p power of Minkowski distance. See \code{\link{dist}} for
#'     details.
#'
#' @param stand logical indicating whether the variables should be
#'     standardized for the distance calculation (default: TRUE).
#'
#' @param method the agglomeration method to be used. See
#'     \code{\link{hclust}} for details.
#' 
#' @param format indicates in which format community composition is
#'     passed. Can be either \code{default} or \code{matrix}.
#' 
#' @param which indicates which results to return
#'     (\code{default="FD"}). The vector may contain any combination
#'     of \code{"FD"},\code{"stdFD"}, \code{"maxFD"},
#'     \code{"distances"}, \code{"tree"},
#'     \code{"branches"}. \code{"FD"} indicates functional diversity
#'     (i.e. total branch length of subtree), \code{"stdFD"}
#'     standardized functional diversity (branch length divided by
#'     total tree branch length) and \code{"maxFD"} total tree branch
#'     length. \code{"distances"} also returns the calculated distance
#'     matrix, \code{"tree"} the calculated trait tree, and
#'     \code{"branches"} a list of branch indiced within the tree that
#'     connect a tip to the root.
#' 
#' @param shrink.tree logical indicating whether species not present
#'     in any community should be removed from the trait
#'     tree. defaults to \code{TRUE}.
#' 
#' @return a list containing the items selected by \code{which}.
#' 
#' @seealso \code{\link{calcPD}} 
#' @examples
#' library(pdiv)
#' 
#' data(composition.example)
#' composition.example
#' ##    community species
#' ## 1          1       a
#' ## 2          1       b
#' ## 3          1       c
#' ## ...
#' ## 19         7       c
#' ## 20         8       d
#' 
#' data(traits.example)
#' traits.example
#' ##    species     trait1    trait2     trait3
#' ## 1        a 0.39486189 0.9727755 0.83947128
#' ## 2        b 0.37285983 0.6587653 0.94541045
#' ## 3        c 0.06773421 0.0873821 0.39090723
#' ## ...  
#' ## 9        i 0.95826170 0.1678766 0.92731325
#' ## 10       j 0.29947903 0.3547652 0.36046854
#' 
#' calcFD(composition.example,traits.example)
#' ## $FD
#' ##    1    2    3    4    5    6    7    8
#' ## 1.63 1.25 1.46 2.10 0.00 0.00 0.00 0.00 
#'
#' @references Petchey O \& Gaston K (2006) Functional diversity: back to basics and looking forward.
#' Ecology Letters 9: 741--758.
#'
#' Petchey OL \& Gaston KJ (2002) Functional diversity (FD), species
#' richness and community composition. Ecology Letters 5: 402--411.
#'
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom cluster daisy
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' @importFrom ape as.phylo
#' @export    
calcFD <- function(communities = NULL, traits=NULL, distance="euclidean", p=2, stand=TRUE,
                   method="complete", format="default", which=c("FD"),shrink.tree=TRUE) 
{
    requireNamespace("cluster")
    rr <- list();
    which.options <- c("distances","tree","branches","FD","stdFD","maxFD");
    which <- which.options[ pmatch(which, which.options) ];

    formats <- c("default", "matrix");
    format <- pmatch(format, formats)
    if (is.na(format)) 
        stop("invalid format for community composition")
    if (format == -1) 
        stop("ambiguous format for community composition")

    communities <- as.data.frame(communities);  
    if( format == "matrix" ) {
        rownames( communities ) <- communities[ , 1 ]
        communities <- communities[ , -1 ]
    }
    
    if("matrix" %in% format) 
        communities <- communityFromMatrix( communities );
    
    sp.col.trt <- which( colnames(traits) %in% colnames(communities) )
    if( length(sp.col.trt) == 0 ) 
        stop("Names of columns containing species labels do not match")
    if( length( sp.col.trt ) != 1 ) 
        stop("Name of column containing species labels ambiguous")
    
    sp.colname  <- colnames( traits )[ sp.col.trt ]
    sp.col.com  <- which( colnames(communities) == sp.colname )
    com.col.com <- 3 - sp.col.com

    rownames(traits) <- as.character(traits[,sp.col.trt]);

    ## remove all species from trait tree that are not in a community
    if(shrink.tree) {        
        sp.to.keep <- which(!is.na( 
            match(rownames(traits),
                  as.character(communities[,sp.col.com]) )
            ))
        traits<-traits[sp.to.keep,-sp.col.trt,drop=FALSE];
    } else {
        traits<-traits[,-sp.col.trt,drop=FALSE];
    }

    ## calculate trait tree
    distances <- if(distance == "gower")
                     daisy(traits, metric = "gower", stand = stand)
                 else 
                     dist( if(stand) scale(traits) else traits,  method = distance )
    if("distances" %in% which)
        rr$distances <- distances
    tree <- as.phylo( hclust( distances, method=method ) )
    if("tree" %in% which)
        rr$tree <- tree;
    
    ## walk up tree from all tips, collect edge indices
    ntip <- length(tree$tip.label);
    branches <- lapply(1:ntip, function(x) .walk_to_root(tree,x) ); 

    if("branches" %in% which)
        rr$branches <- branches;

    names(branches) <- tree$tip.label
    branches <- .remove_common(branches);
    FDmax <- sum(tree$edge.length[unique(unlist(branches))], na.rm=TRUE)
    if("maxFD" %in% which)
        rr$maxFD <- FDmax;

    if("FD" %in% which || "stdFD" %in% which) {
        r <- data.frame(com = sort(unique(as.character(communities[,com.col.com]))), FD = NA)
        for(i in 1:nrow(r)) {
            sp.set <- sort(unique(communities[communities[,com.col.com] == r$com[i],sp.col.com]));
            br <- branches;

            if(!all(sp.set %in% names(br))) 
                stop("Species not found in trait tree: ",
                     paste(sp.set[which(!(sp.set %in% names(br)))],
                           collapse=", "));

            br <- br[match(sp.set,names(br))]  
            br <- .remove_common(br)
            r$FD[i] <- sum(tree$edge.length[unique(unlist(br))],na.rm=TRUE)
        }
        if("FD" %in% which) {
            rr$FD <- r$FD;
            names(rr$FD) <- r$com;
        }
        if("stdFD" %in% which) {
            rr$stdFD <- r$FD/FDmax;
            names(rr$stdFD) <- r$com;
        }
    }
    rr;
}
