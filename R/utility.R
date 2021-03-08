##' Sets of form HuA for A subset of T
##'
##' @param h,t head and tail sets
##'
##' @details Returns a list of sets of the form
##' \eqn{H \cup A} where \eqn{A \subseteq T}.
##'
tailPowSet <- function(h,t) {
  ps <- rje::powerSet(t)
  lapply(ps, function(x) sort.int(c(h,x)))
}

##' Assign signs to head-tail sets
##'
##' @param h,t head and tail sets
##'
##' Assign -1 to sets of even size, 1 to sets of odd size.
setSign <- function(h,t) {
  out <- c(1)
  for (i in seq_along(t)) out <- kronecker(out, c(1,-1))

  # (-1)^(length(h)+1)*rje::subsetMatrix(length(t))[,1]
  (-1)^(length(h)+1)*out
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
##' wh_entries(list(1, 3, 1:3))
##' wh_entries(powerSet(1:3))
##' wh_entries(1:3)
##'
wh_entries <- function(x) {
  if (!is.list(x)) x <- list(x)
  if (length(x) == 0) return(integer(0))

  ## get maximum value
  n <- max(unlist(x))

  wgts <- 2^(seq_len(n)-1)

  ## return locations in imset vector
  sapply(x, function(x) sum(wgts[x])) + 1
}
