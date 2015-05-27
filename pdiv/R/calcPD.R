#' Calculate phylogenetic diversity 
#'
#' There are also other routines to calculate phylogenetic
#' diversity. I implemented this function to have one that follows
#' exactly the same procedure as the metric defined by Petchey &
#' Gaston (see \code{\link{calcFD}}). Phylogenetic diversity is
#' calculated as total branch length of the phylogenetic subtree
#' defined by the set of species present in a community.
#' 
#' @param communities data frame containing a column with species
#' names and a column designating the communities (e.g. plots or other
#' groups to compare)
#' @param format indicates in which format community composition is
#' passed. Can be either \code{default} or \code{matrix}
#' @param tree phylogenetic tree of class \code{phylo}
#' @param which indicates which results to return
#' (\code{default="PD"}). The vector may contain any combination of
#' \code{"PD"}, \code{"maxPD"}, and \code{"branches"}. \code{"PD"}
#' indicates phylogenetic diversity (i.e. total branch length of
#' subtree), \code{"maxPD"} total tree branch length, and
#' \code{"branches"} a list of branch indiced within the tree that
#' connect a tip to the root.
#' 
#' @return a list containing the items selected by \code{which}.
#' @seealso \code{\link{calcFD}} 
#' @examples
#'
#' library(pdiv)
#' data(composition.example)
#' data(phylotree.example)
#' plot(phylotree.example)
#' calcPD(composition.example,phylotree.example)
#' ## $PD
#' ##   com        PD
#' ## 1   1 1315.1436
#' ## 2   2  254.8177
#' ## 3   3  878.9794
#' ## 4   4 1520.9276
#' ## 5   5    0.0000
#' ## 6   6    0.0000
#' ## 7   7    0.0000
#' ## 8   8    0.0000
#'
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' @export 
calcPD <- function(communities = NULL, tree=NULL, which=c("PD"), format="default") 
{
    rr<-list();

    which.options <- c("branches","PD","maxPD");
    which <- which.options[ pmatch( which, which.options ) ]

    formats <- c("default", "matrix");
    format <- pmatch(format, formats)
    if (is.na(format)) 
        stop("invalid format for community composition");
    if (format == -1) 
        stop("ambiguous format for community composition");      

    communities <- as.data.frame(communities);  
    if(format %in% "matrix") 
        communities <- communityFromMatrix(communities);

    all.sp <- sort(unique(as.character(communities[,2])));
    if(!(all(all.sp %in% tree$tip.label))) 
        stop("Species",paste(all.sp[!(all.sp %in% tree$tip.label)],collapse=", "),"missing in tree");

    ## walk up tree from all tips, collect edge indices
    ntip <- length(tree$tip.label);
    branches <- lapply(1:ntip,function(x) .walk_to_root(tree,x)); 
    names(branches) <- tree$tip.label
    branches <- .remove_common(branches);
    PDmax <- sum(tree$edge.length[unique(unlist(branches))],na.rm=T)

    if("branches" %in% which)
        rr$branches <- branches;
    
    if("maxPD" %in% which)
        rr$maxPD <- PDmax;

    if("PD" %in% which) {  
        r <- data.frame(com=sort(unique(communities[,1])),PD=NA)
        for(i in 1:nrow(r)) {
            sp.set <- sort(unique(communities[communities[,1]==r$com[i],2]));
            br <- branches[match(sp.set,names(branches))];
            br <- .remove_common(br);
            r$PD[i] <- sum(tree$edge.length[unique(unlist(br))],na.rm=T);
        }
        rr$PD <- r;
    }
    rr;
}
