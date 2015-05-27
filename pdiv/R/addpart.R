#' Additive Partitioning
#'
#' This function computes Loreau and Hectors additive partitioning,
#' given a data set that specifies the composition of plant
#' communities on a species basis.
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
#' @param d Data frame describing the composition of plant
#' communities, on a per-species basis
#' @param depmix A 'model formula' describing the structure of the
#' data, in the \code{form y ~ mixture/species + plot}
#' @param name.CE Name of column holding the calculated complemetarity
#' effect. Defaults to \code{CE.y}, where \code{y} is replaced by the
#' name of the dependent variable (usually the biomass measure).
#' @param name.SE Name of column holding the calculated selection
#' effect. Defaults to \code{SE.y}, where \code{y} is replaced by the
#' name of the dependent variable specified.
#' @param excl.zeroes Exclude species with unknown or zero biomass
#' from calculations.
#' @return A data frame containing columns for plant species mixture,
#' plot, and the computed complementarity and selection effects.
#' 
#' @examples
#'
#' ## TO BE ADDED
#' 
#' @references Loreau M \& Hector A (2001) Partitioning selection and
#' complementarity in biodiversity experiments. Nature 412, 72-76
#'
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' @export    
addpart <- function(d, depmix, name.CE=NULL,name.SE=NULL, excl.zeroes=FALSE)
{
    ## define population covariance
    covp <- function(x,y) {
        sum((x-mean(x))*(y-mean(y)))/length(x);
    }

    ## check for formula structure
    if(length(depmix)!=3
       || length(depmix[[3]][[2]])!=3
       || length(depmix[[3]][[3]])!=1) {
        stop("formula <",
             deparse(depmix),
             "> does not conform to required structure, ",
             "which is y ~ mixture/species + plot",sep="");
    }

    yname   <- as.list(depmix)[[2]];
    mixname <- as.character(as.list(depmix)[[3]][[2]][[2]]);
    spname  <- as.character(as.list(depmix)[[3]][[2]][[3]]);
    unitname<- as.character(as.list(depmix)[[3]][[3]]);
    if(is.null(name.CE))
        name.CE <- paste("CE",yname,sep=".");
    if(is.null(name.SE))
        name.SE <- paste("SE",yname,sep=".");

    d0      <- model.frame(depmix,data=d);
    col.y   <- which(names(d0)==yname); 
    col.mix <- which(names(d0)==mixname);
    col.sp  <- which(names(d0)==spname);
    col.unit<- which(names(d0)==unitname);

    d0[[col.mix]] <- as.character(d0[[col.mix]]);
    d0[[col.sp]]  <- as.character(d0[[col.sp]]);
    d0[[col.unit]]<- as.character(d0[[col.unit]]);

    ## identify monocultures
    d0mono <- sapply(d0[[ col.mix ]],
                     function(mixref) {
                         0 == sum( sapply(d0[[ col.mix ]],
                             function(mix) {
                                 mixref != mix & grepl(mix,mixref,fixed=TRUE) }
                                          ))
                     })
    
    ## prepare empty data frame for CE, SE
    sepchar <- '\01'

    comb <- unique(paste(d0[[col.mix]],d0[[col.unit]],sep = sepchar)); 
    r <- data.frame(sapply(strsplit(comb,sepchar),function(x) x[1]),
                    sapply(strsplit(comb,sepchar),function(x) x[2]));

    names(r) <- c(mixname,unitname);

    ## now loop over all mixtures and determine CE and SE
    r[[name.CE]]<-NA;
    r[[name.SE]]<-NA;
    for(i in seq(along=rownames(r))) {
        ## extract all species contained
        sp <- sort(unique(as.character(d0[[col.sp]][d0[[col.mix]]==r[i,1]])));
        S <- length(sp);
        
        if(S>1) {
            ## collect biomass in mixtures
            m <- sapply(sp,
                        function(sp) {
                            mean(d0[  d0[[col.mix]]==r[i,1] 
                                    & d0[[col.unit]]==r[i,2] 
                                    & d0[[col.sp]]==sp, col.y])
                        });

            ## collect biomass in monocultures
            m0 <- sapply(sp,
                         function(sp) {
                             mean(d0[d0mono
                                     & grepl(sp,d0[[col.mix]],
                                             fixed=TRUE),
                                     col.y],
                                  na.rm=TRUE)
                         })
            
            if(excl.zeroes) {  
                idx <- (is.na(m0) | m0==0);          
                S<-S-sum(idx);
                m<-m[!idx];
                m0<-m0[!idx];
            }
            
            if(any(m0==0) || any(is.na(m0))) {
                warning("Could not calculate CE/SE because at least ",
                        "one monoculture biomass is zero, NA, or missing");
            } else {
                deltaRY <- m/m0 - 1/S;
                if(length(deltaRY)>1) {
                    r[[i,name.CE]] <- S*mean(deltaRY)*mean(m0);
                    r[[i,name.SE]] <- S*covp(deltaRY,m0);
                }
            }
        } 
    }
    return(r);
}
