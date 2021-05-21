##' Compute subsets by integer position in vector
##'
##' @param x vector of integers to compare
##' @param simplify if \code{x} has length one, just return a vector
##'
##' @export
posToSubset <- function(x, simplify=FALSE) {
  if (simplify && length(x) == 1) {
    return(which(as.numeric(intToBits(x-1)) > 0))
  }
  lapply(x, function(z) which(as.numeric(intToBits(z-1)) > 0))
}

##' Compare subsets by integer position in vector
##'
##' @param x,y values to compare
##'
##' @export
intSubset <- function(x, y) {
  all(as.numeric(intToBits(y-1)) - as.numeric(intToBits(x-1))) >= 0
}
