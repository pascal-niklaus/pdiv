#' Additive and Tripartite Partitioning
#'
#' Given a data set that specifies the composition of plant
#' communities on a per species basis, these functions compute Loreau and
#' Hector's additive partitioning, or Fox' tripartite partitioning.
#'
#' The input data frame needs to contain columns defining the units
#' for which complementarity and selection effects are to be
#' calculated (typically field plots), the community composition of
#' these plots, the species contained in these plots, plus their
#' contribution to the community-level metric that is analysed
#' (typically biomass). These columns are specified using a model
#' formula of the form \code{form y ~ mixture/species + plot}. Note
#' that plot and mixture can be identical, in which case one has to
#' write \code{form y ~ mixture/species + mixture}.
#'
#' Species composition is a character string constructed so that the
#' presence of individual species can be tested by searching for the
#' occurence of the species code in the composition string. For
#' example, species \code{A} will be found in composition \code{ABCD}
#' but not in \code{EFGH}. The composition codes can contain
#' additional characters as long as the above rule is met,
#' e.g. \code{A|B|C|D} would work as well.
#'
#' Note: Monocultures are identified by testing if their composition
#' codes do not contain any of the other composition codes. It
#' therefore is crucial to include the monocultures in the data set,
#' even if their biomass is zero or not known. Otherwise, mixtures may
#' be misidentified as monocultures (because there are no
#' corresponding lines for their component monocultures). Example:
#' composition 'AB' will be treated as monoculture if no composition
#' 'A' or 'B' are present in the data set.
#' 
#' Monoculture biomass (the reference for the underlying relative
#' yield calculation) is calculated as mean of all the monoculture
#' plots, excluding \code{NA}s.
#'
#' Relative yields, on which the additive partitioning scheme bases,
#' cannot be calculated for communities that contain species with zero
#' or unknown monoculture biomass. For such plots, the complementarity
#' and selection effects will be set to \code{NA}, unless
#' \code{excl.zeroes=TRUE} is specified. Then, these species will be
#' removed from the species set before complementarity and selection
#' effects are calculated. The underlying relative yield data are
#' still calculated using the original species numbers; only then
#' are the calculated CE plus SE equal to the net effect calculated
#' as mixture yield minus average monoculture yield (incl. the failing
#' monoculture with zero yield).
#'
#' These partitioning schemes can be performed by groups
#' (e.g. different experimental treatments or dates). These groups are
#' specified as right-hand side only formula. If the code runs on a
#' platform other than Windows and the library \code{parallel} is
#' installed, the processing of the different groups is parallelized
#' on different processor cores.
#'
#' @param data Data frame describing the composition of plant
#'     communities, on a per-species basis
#' 
#' @param depmix A 'model formula' describing the structure of the
#'     data, in the \code{form y ~ mixture/species + plot}
#'
#' @param name.CE Name of column holding the calculated complemetarity
#'     effect. Defaults to \code{CE.y}, where \code{y} is replaced by
#'     the name of the dependent variable (usually the biomass
#'     measure).
#'
#' @param name.SE Name of column holding the calculated selection
#'     effect. Defaults to \code{SE.y}, where \code{y} is replaced by
#'     the name of the dependent variable specified.
#'
#' @param name.TIC Name of column holding the calculated
#'     trait-independent complemetarity effect. Defaults to
#'     \code{TIC.y}, where \code{y} is replaced by the name of the
#'     dependent variable (usually the biomass measure). Numerically
#'     this term is identical to the complementarity effect of the
#'     additive partitioning scheme.
#'
#' @param name.DOM Name of column holding the calculated selection
#'     dominance effect. Defaults to \code{DOM.y}, where \code{y} is
#'     replaced by the name of the dependent variable specified.
#'
#' @param name.TDC Name of column holding the calculated
#'     trait-dependent complementarity effect. Defaults to
#'     \code{TDC.y}, where \code{y} is replaced by the name of the
#'     dependent variable specified.
#'
#' @param excl.zeroes logical specifying whether species with unknown
#'     or zero biomass are excluded from calculations (default \code{FALSE}).
#'     Data that are NA are also excluded.
#'
#' @param groups The calculations can be performed by groups specified
#'     with a right hand side only model formula.
#'
#' @return A data frame containing columns for plant species mixture,
#'     plot, and the computed complementarity and selection effects.
#'
#' @examples
#'
#' data(divpart)
#' # additive partitioning
#' addpart(m ~ comp/sp+plot, data = divpart, groups = ~group)
#' # tripartite partitioning
#' tripart(m ~ comp/sp+plot, data = divpart, groups = ~group)
#'
#' @references Loreau M \& Hector A (2001) Partitioning selection and
#'     complementarity in biodiversity experiments. Nature 412, 72-76
#'
#' @references Fox (2005) Interpreting the selection effect of
#'     biodiversity on ecosystem functioning. Ecology Letters 8,
#'     846--856
#'
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#'
#' @importFrom stats model.frame na.pass
#' @rdname divpartfun
#' @export
addpart <- function(depmix,
                    data,
                    name.CE = NULL,
                    name.SE = NULL,
                    excl.zeroes = FALSE,
                    groups = NULL)
{
    .divpart(depmix = depmix,
             d = data,
             name.CE = name.CE,
             name.SE = name.SE,
             excl.zeroes = excl.zeroes,
             groups = groups,
             scheme = "ap")
}

#' @rdname divpartfun
#' @export
tripart <- function(depmix,
                    data,
                    name.TIC = NULL,
                    name.TDC = NULL,
                    name.DOM = NULL,
                    excl.zeroes = FALSE,
                    groups = NULL)
{
    .divpart(depmix = depmix,
             d = data,
             name.TIC = name.TIC,
             name.TDC = name.TDC,
             name.DOM = name.DOM,
             excl.zeroes = excl.zeroes,
             groups = groups,
             scheme = "tp")
}

## private function called by addpart and tripart
## (note that depmix and groups are evaluated in the frame of the caller
## of addpart and tripart, i.e the parent's parent frame)
.divpart <- function(depmix,
                     d,
                     name.CE = NULL,
                     name.SE = NULL,
                     name.TIC = NULL,
                     name.TDC = NULL,
                     name.DOM = NULL,
                     excl.zeroes, groups,
                     scheme)
{
    ## evaluate formula in caller context of addpart and tripart
    if(!is.null(groups)) {
        if(class(groups) != "formula" || length(groups) != 2)
            stop("argument 'groups' (",deparse(groups),
                 ") should be a right-hand side formula")
        d.grouping <- eval(model.frame(groups, d, na.action = na.pass), parent.frame(2))
    }

    d <- eval(model.frame(depmix, d, na.action=na.pass), parent.frame(2))
    
    ## create list with data groups (dlist)
    if(is.null(groups))
        dlist <- list(all = d)
    else
        dlist <- split(
            cbind(d.grouping,d),
            f = .columnsToNum(d.grouping))

    ## check for formula structure
    if(class(depmix) != "formula" ||
       length(depmix) != 3 ||
       length(depmix[[3]][[2]]) != 3 ||
       length(depmix[[3]][[3]]) != 1)
        stop("formula <",
             deparse(depmix),
             "> does not conform to required structure, ",
             "which is y ~ mixture/species + plot",sep="")
    yname   <- as.character(as.list(depmix)[[2]])
    mixname <- as.character(as.list(depmix)[[3]][[2]][[2]])
    spname  <- as.character(as.list(depmix)[[3]][[2]][[3]])
    unitname<- as.character(as.list(depmix)[[3]][[3]])

    ## create result columns depending on partitioning scheme
    if(scheme == "ap") {
        if(is.null(name.CE))
            name.CE <- paste("CE", yname, sep=".");
        if(is.null(name.SE))
            name.SE <- paste("SE", yname, sep=".");
    } else if(scheme == "tp") {
        if(is.null(name.TIC))
            name.TIC <- paste("TIC", yname, sep=".")
        if(is.null(name.TDC))
            name.TDC <- paste("TDC", yname, sep=".")
        if(is.null(name.DOM))
            name.DOM <- paste("DOM", yname, sep=".")
    } else
        stop("Unknown partitioning scheme")

    col.y   <- which(names(dlist[[1]]) == yname)
    col.mix <- which(names(dlist[[1]]) == mixname)
    col.sp  <- which(names(dlist[[1]]) == spname)
    col.unit<- which(names(dlist[[1]]) == unitname)

    ## process data groups
    r <-
        .mylapply(
            dlist,
            function(d) {
                if(!is.null(groups)) {
                    grpcodes <- model.frame(groups,d)[1,,drop=FALSE]
                    rownames(grpcodes) <- NULL
                }
                d[[col.mix]] <- as.character(d[[col.mix]])
                d[[col.sp]]  <- as.character(d[[col.sp]])
                d[[col.unit]]<- as.character(d[[col.unit]])
                
                ## identify monocultures and perform sanity checks
                ## monocultures are the compositions that do not contain
                ## any other composition code
                comps <- sort(unique(d[[col.mix]]))

                monos <-
                    sapply(
                        comps,
                        function(mixref)
                            0 == sum(sapply(comps,
                                            function(mix)
                                                (mixref != mix) &&
                                                    grepl(mix, mixref, fixed=TRUE)))
                    )

                rich <-
                    sapply(
                        comps,
                        function(mix)
                            length(unique(d[[col.sp]][ d[[col.mix]] == mix ] ))
                    )

                if(any(rich[monos] != 1)) 
                    stop("Monoculture(s) ",
                         paste(paste("'",comps[monos & rich>1], "'", sep="", collapse=", ")),
                         " contain(s) multiple species.")
                
                if(any(rich[!monos] == 1))
                    warning("Mixture(s) ",
                            paste(paste("'",comps[!monos & rich==1], "'", sep="", collapse=", ")),
                            " contain(s) only one species.")
                
                dmono <- d[[ col.mix ]] %in% comps[monos]

                ## prepare empty data frame for CE, SE
                r <- unique(d[,c(col.mix,col.unit)])                
                names(r) <- c(mixname,unitname)

                if(!is.null(groups))
                    r <- data.frame(grpcodes, r)
                if(scheme == "ap") {
                    r[[name.CE]] <- NA
                    r[[name.SE]] <- NA
                } else if(scheme == "tp") {
                    r[[name.TIC]] <- NA
                    r[[name.TDC]] <- NA
                    r[[name.DOM]] <- NA
                }

                ## loop over all mixtures and determine CE and SE
                for(i in seq(along = rownames(r))) {
                    ## extract all species contained in mixture
                    sp <- sort(unique(as.character(
                        d[[ col.sp ]][ d[[col.mix]] == r[i, mixname] ])))
                    S <- length(sp)
                    if(S > 1) {
                        ## collect biomass in mixtures
                        m <- sapply(sp,
                                    function(sp) {
                                        mean(d[ d[[col.mix]] == r[i, mixname] &
                                                d[[col.unit]] == r[i, unitname] &
                                                d[[col.sp]] == sp, col.y ] )
                                    })
                        ## collect biomass in monocultures
                        m0 <- sapply(sp,
                                     function(sp) {
                                         mean(d[ dmono &
                                                 grepl(sp, d[[col.mix]],
                                                       fixed = TRUE),
                                                col.y ],
                                              na.rm = TRUE)
                                     })

                        idx <- excl.zeroes & ( is.na(m0) | ( m0 == 0) )                        
                        Sreduced <- S - sum( idx )
                        m <- m[ ! idx ]
                        m0 <- m0[ ! idx ]
                        
                        if( any(m0 == 0) || any(is.na(m0)) ) {
                            warning("Could not calculate partitioning for ",r[i, mixname],
                                    " because at least ",
                                    "one monoculture biomass is zero, NA, or missing");
                        } else {
                            RY <- m / m0
                            RYE <- 1 / S
                            deltaRY <- RY - RYE
                            if(length(deltaRY ) > 1) {
                                if(scheme == "ap") {
                                    r[[i, name.CE]] <- Sreduced * mean(deltaRY) * mean(m0)
                                    r[[i, name.SE]] <- Sreduced * .covp(deltaRY, m0 )
                                } else if(scheme == "tp") {
                                    RYT <- sum( RY )
                                    r[[i, name.TIC]] <- Sreduced * mean(deltaRY) * mean(m0)
                                    r[[i, name.DOM]] <- Sreduced * .covp(RY / RYT - RYE, m0)
                                    r[[i, name.TDC]] <- Sreduced * .covp(RY - RY / RYT, m0)
                                }
                            }
                        }
                    }
                }
                r
            }
        )
    ## combine all groups
    r <- r[order(names(r))]
    r <- do.call("rbind",r)
    rownames(r) <- NULL
    r
}
