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

##' Vector of sizes for power set
##'
##' @param n size of set
##' @param m maximum size to consider
##'
##' @details Uses Kronecker sum to obtain sizes of subsets in a conventionally
##' ordered power set.
##'
setSize <- function(n, m=n) {
  out <- 0
  for (i in seq_len(n)) out <- kronecker(out, c(0,1), FUN="+")

  out <- out[out <= m]

  out
}


##' @describeIn elem_imset
##'
##' @param ci object of class \code{ci}
##' @param imset semi-elementary imset
##'
##' @export
sei2ci <- function(imset) {

  tmp <- posToSubset(which(imset != 0))

  if (length(tmp) == 0) return(as.ci(list(integer(0),integer(0))))
  else if (length(tmp) != 4) stop("Not a semi-elementary imset")

  out <- list()

  A <- setdiff(tmp[[2]], tmp[[1]])
  B <- setdiff(tmp[[3]], tmp[[1]])
  C <- tmp[[1]]

  out <- as.ci(list(A,B,C))
  if (!setequal(tmp[[4]], c(A,B,C))) stop("Not a semi-elementary imset")

  out
}

##' @describeIn elem_imset convert conditional independence to semi-elementary imset
##' @export
ci2sei <- function(ci, n=max(vars)) {
  vars <- unlist(ci)

  elem_imset(ci[[1]], ci[[2]], ci[[3]], n=n)
}
