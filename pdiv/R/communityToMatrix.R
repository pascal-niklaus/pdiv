#' Interconvert formats of community composition data
#'
#' Community composition can be represented in different
#' formats. These function interconvert matrix-type representations
#' into data frames with communitiy and species in parallel columns.
#' 
#' @param cc data specifying community composition. For
#' \code{communityFromMatrix}, this is a matrix specifying the
#' presence of species. Absence of species is indicated by zero or NA,
#' any other value is interpreted as species being present.  For
#' \code{communityToMatrix}, this is a data frame listing the species
#' present by means of two separate columns for community and species,
#' plus an optional column with abundances.
#' @param column.with.community column containing the community code,
#' specified as column number. Defaults to \code{1}. Can be set to
#' \code{NA} if the community code is stored as row name.
#' @param column.with.species column containing the species codes,
#' either specified as column number of name. Defaults to \code{2}.
#' @param column.with.abundances column with species abundances, if
#' these are provided (default \code{NULL}).
#' @param columns.with.species Columns of data matrix containing
#' species presence data. Specified in numeric format. Defaults to all
#' columns except the one containing the community code.
#' @param community.in.column logical. Definies whether the community
#' code should be stored as separate column of the resulting data
#' frame (default, \code{TRUE}), or as row names (\code{FALSE}).
#' @param community.column.name name of the data column containing the
#' community codes, in case \code{community.in.column=TRUE}. Defaults
#' to \code{"community"}.
#' @param species.column.name name of the data column containing the
#' species codes. Defaults to \code{"species"}.
#' @param abundance.column.name name of column containing abundance
#' data, or NULL if none
#' @param species.missing value in matrix indicating missing
#' species. Defaults to \code{0}, but could for example also be set to
#' \code{NA}.
#' @param keep.abundances logical determining whether abundances are converted
#' to presence/absence data (default \code{FALSE})
#' @return A data frame containing community composition in the respective format
#' @examples
#' library(pdiv)
#' data(composition.example)
#' ##    community species
#' ##            1       a
#' ##            1       b
#' ##            1       c
#' ##            1       f
#' ##            2       f
#' ##          ...     ...
#' ##            7       c
#' ##            8       d
#' cmatrix <- communityToMatrix(composition.example)
#' cmatrix
#' ##   community a b c d e f h j
#' ## 1         1 1 1 1 0 0 1 0 0
#' ## 2         2 0 0 0 0 0 1 1 1
#' ## 3         3 1 0 0 0 1 0 0 1
#' ## 4         4 1 1 1 1 1 0 1 0
#' ## 5         5 1 0 0 0 0 0 0 0
#' ## 6         6 0 1 0 0 0 0 0 0
#' ## 7         7 0 0 1 0 0 0 0 0
#' ## 8         8 0 0 0 1 0 0 0 0
#' communityFromMatrix(cmatrix)
#' ##    community species
#' ## 1          1       a
#' ## 2          1       b
#' ## 3          1       c
#' ## ...      ...     ...
#' ## 19         7       c
#' ## 20         8       d
#'
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' @rdname communityToMatrix
#' @export    
communityToMatrix <- function(
    cc,
    column.with.community = 1,
    column.with.species = 2,
    column.with.abundances = NULL,
    keep.abundances = TRUE,
    community.in.column = TRUE,
    community.column.name = "community",
    species.missing = 0)
{
    abu <- if(is.null(column.with.abundances)) {
        rep(1,nrow(cc))
    } else {
        cc[,column.with.abundances]
    }
    if(!keep.abundances)
        abu <- sign(abu)
    com <- sort(unique(cc[,column.with.community]))
    sp  <- sort(unique(cc[,column.with.species]))
    x <- as.data.frame(
        matrix(data=species.missing,
               nrow=length(com),
               ncol=length(sp)))
    colnames(x) <- sp
    rownames(x) <- com
    for(i in 1:nrow(x)) {
        cidx <- cc$com == com[i]
        x[i,] <- 
            sapply(match(sp,cc$species[cidx]),
                   function(x) if(is.na(x)) species.missing else abu[cidx][x])
    }
    if(community.in.column) {
        x[[community.column.name]] <- row.names(x)
        x <- x[,c(ncol(x),1:(ncol(x)-1))]
    }
    x
}

#' @rdname communityToMatrix
#' @export    
communityFromMatrix <- function(
    cc,
    column.with.community = 1,
    columns.with.species = setdiff(1:ncol(cc),column.with.community),
    community.column.name =
        if(is.na(column.with.community))
            "community"
        else
            colnames(cc)[column.with.community],
    species.column.name = "species",
    abundance.column.name = NULL)
{
    present <- function(x) {
        !(is.na(x) | x==0)
    }
    x <- data.frame(community = rep( if(is.na(column.with.community))
                                         rownames(cc)
                                     else
                                         cc[,column.with.community],
                        each=length(columns.with.species)), 
                    species = colnames(cc)[columns.with.species])
    x$abundance <- c(t(as.matrix(cc[,columns.with.species])))
    x <- x[ present(x$abundance), ]
    if(is.null(abundance.column.name)) {
        x <- x[,-3]
    }
    x <- x[order(x$community,x$species),]
    rownames(x) <- 1:nrow(x)
    colnames(x) <- c(community.column.name,
                     species.column.name,
                     abundance.column.name)
    x
}
