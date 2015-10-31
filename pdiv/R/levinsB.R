#' Levins' B
#'
#' Calculate niche breadth as defined by Levins (1968).
#'
#' @param x Measure of species abundance in niches. These need not
#' be standardized to unity sum.
#' @return Niche breadth as defined by Levins.
#' 
#' @references Levins (1968) Evolution in changing
#'     environments. Princeton University Press, Princeton, New
#'     Jersey, USA.
#' 
#' @examples
#' library(pdiv)
#' x <- c(0,1,2,3,2,0)
#' levinsB(x)
#' 
#' @author Pascal Niklaus \email{pascal.niklaus@@ieu.uzh.ch}
#' 
#' @export  
levinsB <- function(x) 
{ 
  sum(x)^2/(sum(x^2));
}
