##' Move between subsets and integer position in vector
##'
##' @param x vector of integers to compare
##' @param y list of subsets to consider
##' @param simplify if \code{x} has length one, just return a vector
##' @name position_subset
NULL

##' @describeIn position_subset Obtain subset from position
##' @export
posToSubset <- function(x, simplify=FALSE) {
  if (simplify && length(x) == 1) {
    return(which(as.numeric(intToBits(x-1)) > 0))
  }
  lapply(x, function(z) which(as.numeric(intToBits(z-1)) > 0))
}

##' @describeIn position_subset Integer position from subset
##' @export
subsetToPos <- function(y) {
  sapply(y, function(k) sum(2^(k-1))+1)
}

##' Get entries in vector of subsets
##'
##' Gives location where set would appear in imset under
##' lexicographic ordering.
##'
##' @param x list of positive integer subsets
##'
##' Given a subset of 1,2,3,... returns the location
##' where set would appear in lexicographic order.
##'
##' @examples
##' imset:::wh_entries(list(1, 3, 1:3))
##' imset:::wh_entries(powerSet(1:3))
##' imset:::wh_entries(1:3)
##'
wh_entries <- function(x) {
  .Deprecated(subsetToPos, msg="Function deprecated, use subsetToPos() instead")
  if (!is.list(x)) x <- list(x)
  if (length(x) == 0) return(integer(0))

  ## get maximum value
  n <- max(unlist(x))

  wgts <- 2^(seq_len(n)-1)

  ## return locations in imset vector
  sapply(x, function(x) sum(wgts[x])) + 1
}


##' Compare subsets by integer position in vector
##'
##' @param x,y values to compare
##'
##' @export
intSubset <- function(x, y) {
  all(as.numeric(intToBits(y-1)) - as.numeric(intToBits(x-1))) >= 0
}
