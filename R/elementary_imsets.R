##' Get a (semi-)elementary imset
##'
##' Returns the semi-elementary imset \eqn{u_{\langle A, B | C \rangle}}{u<A, B | C>}
##' as defined in Studeny (2006).
##'
##' @param A,B,C disjoint subsets of 1,..,n
##' @param n number of variables involved
##'
##' @details Returns an imset, a vector of length \eqn{2^n}
##' with all but four entries zero.
##'
##' @export
elem_imset <- function(A, B, C=integer(0), n=max(c(A,B,C))) {

  if (length(intersect(A,B)) > 0) stop("sets A and B should not intersect")
  if (length(intersect(c(A,B), C)) > 0) {
    A <- setdiff(A, C)
    B <- setdiff(B, C)
  }
  if (length(A) == 0 || length(B) == 0) return(zero_imset(n))

  ## compute entries in vector
  wts <- 2^(seq_len(n)-1)
  wtC <- sum(wts[C])
  wtA <- sum(wts[A])
  wtB <- sum(wts[B])

  ## construct vector
  out <- rep(0,2^n)
  out[c(wtC, wtC+wtA, wtC+wtB, wtC+wtA+wtB)+1] = c(1,-1,-1,1)

  as.imset(out)
}

##' Deprecated function for (semi-)elementary imsets
##'
##' @param A,B,C disjoint subsets of \eqn{1,...,n}
##' @param n number of variables involved
##'
##' @details Returns an imset, a vector of length \eqn{2^n}
##' with all but four entries zero.
##'
##' Function is now deprecated in favour of \code{elem_imset}.
##'
##' @name elemImset-deprecated
##' @export
elemImset <- function(A, B, C=integer(0), n=max(c(A,B,C))) {
  .Deprecated("elem_imset")

  if (length(intersect(A,B)) > 0) stop("sets A and B should not intersect")
  if (length(intersect(c(A,B), C)) > 0) {
    A <- setdiff(A, C)
    B <- setdiff(B, C)
  }
  if (length(A) == 0 || length(B) == 0) return(zero_imset(n))

  ## compute entries in vector
  wts <- 2^(seq_len(n)-1)
  wtC <- sum(wts[C])
  wtA <- sum(wts[A])
  wtB <- sum(wts[B])

  ## construct vector
  out <- rep(0,2^n)
  out[c(wtC, wtC+wtA, wtC+wtB, wtC+wtA+wtB)+1] = c(1,-1,-1,1)

  as.imset(out)
}

# gr <- graphCr("1->3<->2->4<->1", format = "ADMG")
# gr2 <- graphCr("1->3<-2->4<-1", format = "ADMG")
#
# all.equal(imset(gr), imset(gr2))
#
# imset(graphCr("1,2,3,4", format = "ADMG"))

##' @describeIn elemImset-deprecated Deprecated version of elemToIndep, use \code{sei2ci}
##' @export
elemToIndep <- function(imset) {
  .Deprecated("sei2ci")
  if (all(imset == 0)) return(as.ci(list(integer(0),integer(0))))
  if (sum(imset != 0) != 4) stop("Not an elementary imset")

  ## isolate sets

  wh_cond <- as.integer(intToBits(which(imset > 0)[1]-1))
  wh_indep1 <- as.integer(intToBits(which(imset < 0)[1]-1))
  wh_indep2 <- as.integer(intToBits(which(imset < 0)[2]-1))

  if (any(wh_indep1 - wh_cond < 0) ||any(wh_indep2 - wh_cond < 0)) stop("Not an elementary imset")

  ## construct CI object
  out <- list()
  out[[1]] <- which(wh_indep1 - wh_cond > 0)
  out[[2]] <- which(wh_indep2 - wh_cond > 0)
  out[[3]] <- which(wh_cond > 0)

  class(out) = "ci"
  out
}
