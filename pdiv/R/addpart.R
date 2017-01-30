#' Additive and Tripartite Partitioning
#'
#' Given a data set that specifies the composition of plant
#' communities on a species basis, these functions compute Loreau and
#' Hectors additive partitioning, or Fox' tripartite partitioning.
#'
#' The input data frame needs to contain columns defining the units
#' for which complementarity and selection effects are to be
#' calculated (typically field plots), the community composition of
#' these plots, the species contained in these plots, plus their
#' contribution to the community-level metric that is to be analysed
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
#' Monoculture biomass (the reference for the underlying relative
#' yield calculation) is calculated as mean of all the monoculture
#' plots. 
#' 
#' Relative yields, on which the additive partitioning scheme bases,
#' cannot be calculated for communities that contain species with zero
#' or unknown monoculture biomass. For such plots, the complementarity
#' and selection effects will be set to NA, unless
#' \code{excl.zeroes=TRUE} is specified. Then, these species will be
#' removed from the species set before complementarity and selection
#' effects are calculated.
#'
#' These partitioning schemes can be performed by groups
#' (e.g. different experimental treatments or dates). These groups are
#' specified as right-hand side only formula. The grouping variables
#' must reside in the data frame. If the code runs on a platform other
#' than Windows and the library \code{parallel} is installed, the
#' processing of the different groups is parallelized to different
#' processor cores.
#' 
#' @param data Data frame describing the composition of plant
#'     communities, on a per-species basis
#' @param depmix A 'model formula' describing the structure of the
#'     data, in the \code{form y ~ mixture/species + plot}
#' @param name.CE Name of column holding the calculated complemetarity
#'     effect. Defaults to \code{CE.y}, where \code{y} is replaced by
#'     the name of the dependent variable (usually the biomass
#'     measure).
#' @param name.SE Name of column holding the calculated selection
#'     effect. Defaults to \code{SE.y}, where \code{y} is replaced by
#'     the name of the dependent variable specified.
#' @param name.TIC Name of column holding the calculated
#'     trait-independent complemetarity effect. Defaults to
#'     \code{TIC.y}, where \code{y} is replaced by the name of the
#'     dependent variable (usually the biomass measure). Numerically
#'     this term is identical to the complementarity effect of the
#'     additive partitioning scheme.
#' @param name.DOM Name of column holding the calculated selection
#'     dominance effect. Defaults to \code{DOM.y}, where \code{y} is
#'     replaced by the name of the dependent variable specified.
#' @param name.TDC Name of column holding the calculated
#'     trait-dependent complementarity effect. Defaults to
#'     \code{TDC.y}, where \code{y} is replaced by the name of the
#'     dependent variable specified.
#' @param excl.zeroes Exclude species with unknown or zero biomass
#'     from calculations.
#' @param groups The calculations can be performed by groups specified
#'     with a right hand side only model formula. The grouping
#'     variables must be in the data frame.
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
#' complementarity in biodiversity experiments. Nature 412, 72-76
#' @references Fox (2005) Interpreting the selection effect of biodiversity on
#' ecosystem functioning. Ecology Letters 8, 846--856
#'
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#'
#' @importFrom stats model.frame
#' @rdname divpartfun
#' @export
addpart <- function(depmix,
                    data,
                    name.CE=NULL,name.SE=NULL,
                    excl.zeroes=FALSE, groups=NULL)
{
    .divpart(depmix=depmix,
             d=data,
             name.CE=name.CE,
             name.SE=name.SE,
             excl.zeroes=excl.zeroes,
             groups=groups,
             scheme="ap")
}

#' @rdname divpartfun
#' @export
tripart <- function(depmix,
                    data,
                    name.TIC=NULL,
                    name.TDC=NULL,
                    name.DOM=NULL,
                    excl.zeroes=FALSE,
                    groups=NULL)
{
    .divpart(depmix=depmix,
             d=data,
             name.TIC=name.TIC,
             name.TDC=name.TDC,
             name.DOM=name.DOM,
             excl.zeroes=excl.zeroes,
             groups=groups,
             scheme="tp")
}

## internal function called by addpart and tripart
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
    parallel <-
        (Sys.info()["sysname"] != "Windows") &&
        requireNamespace("parallel")
    myapply <- if(parallel)
        function(x,fun) parallel::mclapply(x,fun,mc.cores = parallel::detectCores())
    else 
        lapply
    
    ## process the groups
    if(is.null(groups))
        dlist <- list(all=d)
    else {
        if(class(groups) != "formula" || length(groups)!=2)
            stop("'groups' should be a right-hand side formula")
        d.grouping <- model.frame(groups,d)
        ## convert groups to numeric to avoid problems with '.' in names
        for(v in names(d.grouping))
            d.grouping[[v]] <- as.numeric(as.factor(d.grouping[[v]]))
        dlist <- split(d, f = d.grouping, sep = ".")
    }
    
    ## check for formula structure
    if(length(depmix)!=3
       || length(depmix[[3]][[2]])!=3
       || length(depmix[[3]][[3]])!=1) {
        stop("formula <",
             deparse(depmix),
             "> does not conform to required structure, ",
             "which is y ~ mixture/species + plot",sep="")
    }    
    yname   <- as.list(depmix)[[2]]
    mixname <- as.character(as.list(depmix)[[3]][[2]][[2]])
    spname  <- as.character(as.list(depmix)[[3]][[2]][[3]])
    unitname<- as.character(as.list(depmix)[[3]][[3]])

    if(scheme=="ap") {
        if(is.null(name.CE))
            name.CE <- paste("CE",yname,sep=".");
        if(is.null(name.SE))
            name.SE <- paste("SE",yname,sep=".");
    } else if(scheme=="tp") {
        if(is.null(name.TIC))
            name.TIC <- paste("TIC",yname,sep=".")
        if(is.null(name.TDC))
            name.TDC <- paste("TDC",yname,sep=".")
        if(is.null(name.DOM))
            name.DOM <- paste("DOM",yname,sep=".")
    } else
        stop("Unknown partitioning scheme")
  
    r <-
        myapply(
            dlist,
            function(d) {
                if(!is.null(groups)) {
                    grpcodes <- model.frame(groups,d)[1,]
                    rownames(grpcodes) <- NULL
                }                
                d0      <- model.frame(depmix,data=d)
                rownames(d0) <- NULL
                col.y   <- which(names(d0) == yname) 
                col.mix <- which(names(d0) == mixname)
                col.sp  <- which(names(d0) == spname)
                col.unit<- which(names(d0) == unitname)

                d0[[col.mix]] <- as.character(d0[[col.mix]]);
                d0[[col.sp]]  <- as.character(d0[[col.sp]]);
                d0[[col.unit]]<- as.character(d0[[col.unit]]);

                ## identify monocultures
                d0mono <-
                    sapply(d0[[ col.mix ]],
                           function(mixref) {
                               0 == sum( sapply(d0[[ col.mix ]],
                                                function(mix) {
                                                    mixref != mix & grepl(mix,mixref,fixed=TRUE) } )) })
                
                ## prepare empty data frame for CE, SE
                sepchar <- '\01'
                comb <- unique(paste(d0[[col.mix]],
                                     d0[[col.unit]],
                                     sep = sepchar))                    
                r <- data.frame(
                    sapply(strsplit(comb,sepchar),function(x) x[1]),
                    sapply(strsplit(comb,sepchar),function(x) x[2]))
                names(r) <- c(mixname,unitname)
                if(!is.null(groups))
                    r <- data.frame(grpcodes,r)
                
                ## now loop over all mixtures and determine CE and SE
                if(scheme=="ap") {
                    r[[name.CE]] <- NA
                    r[[name.SE]] <- NA
                } else if(scheme=="tp") {
                    r[[name.TIC]] <- NA
                    r[[name.TDC]] <- NA
                    r[[name.DOM]] <- NA
                }
                
                for(i in seq(along=rownames(r))) {
                    ## extract all species contained
                    sp <- sort(unique(as.character(d0[[col.sp]][d0[[col.mix]]==r[i,mixname]])))
                    S <- length(sp)                    
                    
                    if(S>1) {
                        ## collect biomass in mixtures
                        m <- sapply(sp,
                                    function(sp) {
                                        mean(d0[  d0[[col.mix]]==r[i,mixname] 
                                                & d0[[col.unit]]==r[i,unitname] 
                                                & d0[[col.sp]]==sp, col.y]) } )

                        ## collect biomass in monocultures
                        m0 <- sapply(sp,
                                     function(sp) {
                                         mean(d0[d0mono
                                                 & grepl(sp,d0[[col.mix]],
                                                         fixed=TRUE),
                                                 col.y],
                                              na.rm=TRUE) } )
                        
                        if(excl.zeroes) {  
                            idx <- (is.na(m0) | m0==0)       
                            S <- S-sum(idx)
                            m <- m[!idx]
                            m0 <- m0[!idx]
                        }
                        
                        if(any(m0==0) || any(is.na(m0))) {
                            warning("Could not calculate partitioning because at least ",
                                    "one monoculture biomass is zero, NA, or missing");
                        } else {
                            RY <- m/m0
                            RYE <- 1/S
                            deltaRY <- RY - RYE
                            if(length(deltaRY)>1) {
                                if(scheme=="ap") {
                                    r[[i,name.CE]] <- S*mean(deltaRY)*mean(m0)
                                    r[[i,name.SE]] <- S*.covp(deltaRY,m0)
                                } else if(scheme=="tp") {
                                    RYT <- sum(RY)
                                    r[[i,name.TIC]] <- S*mean(deltaRY)*mean(m0)
                                    r[[i,name.DOM]] <- S*.covp(RY/RYT-RYE,m0)
                                    r[[i,name.TDC]] <- S*.covp(RY-RY/RYT,m0)
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

## define population covariance
.covp <- function(x,y) {
    sum((x-mean(x))*(y-mean(y)))/length(x);
}
