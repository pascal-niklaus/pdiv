#' Proportional similarity
#'
#' Calculate proportional similarity as defined in Colwell and Futuyma (1971).
#'
#' @param x Measure of species abundance in niches. These need not
#' be standardized to unity sum.
#' @param y Similar to x but for 2nd species
#' @return Proportional similarity index.
#' 
#' @references Colwell and Futuyma (1971) Measurement of niche breadth
#'     and overlap. Ecology 52:567--576.
#' @examples
#' library(pdiv)
#' x <- c(0,0,1,2,1,0)
#' y <- c(0,2,1,0,0,0)
#' prop.sim(x,y)
#' 
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' 
#' @export  
prop.sim <- function(x,y) 
{ 
  x<-x/sum(x);
  y<-y/sum(y);
  1-.5*sum(abs(x-y)); 
}


