# library(MixedGraphs)
# library(ADMGs)

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
##'
wh_entries <- function(x) {
  if (length(x) == 0) return(integer(0))

  ## get maximum value
  n <- max(unlist(x))

  wgts <- 2^(seq_len(n)-1)

  ## return locations in imset vector
  sapply(x, function(x) sum(wgts[x])) + 1
}

# ##' Calculate the ?? imset of an ADMG
# imset <- function(graph) {
#   out <- rep(0,2^(graph$n))
#   out[1] = -1
#   tmp <- headsTails(graph)
#
#   ## get sets and signs
#   sets <- unlist(mapply(tailPowSet, tmp$heads, tmp$tails, SIMPLIFY = FALSE), recursive = FALSE)
#   signs <- unlist(mapply(setSign, tmp$heads, tmp$tails, SIMPLIFY = FALSE))
#
#   ## reorder
#   wts <- 2^(seq(graph$n)-1)
#   wh <- sapply(sets, function(x) sum(wts[x]))
#   out[wh+1] = signs
# #  out[length(out)] = out[length(out)]-1
#
#   return(out)
# }

##' Standard imset (for DAGs only)
##'
##' @param x a graph in ADMGs format or an imset
##' @export
standard_imset <- function(x, ...) {
  UseMethod("standard_imset")
}

setGeneric("standard_imset")

##' Characteristic imset (for ADMGs...)
##'
##' @param x a graph in ADMGs format or an imset
##' @export
char_imset <- function(x, ...) {
  UseMethod("char_imset")
}

setGeneric("char_imset")

##' @method standard_imset mixedgraph
##' @export
standard_imset.mixedgraph <- function(x) {
  if (!is.DAG(x)) stop("This function only currently works for DAGs")

  n <- length(x$vnames)
  out <- rep(0, 2^n)
  out[1] = -1; out[2^n] = 1

  ## get entries associated with parent sets
  pa_set <- sapply(seq_len(n), function(v) {
    sum(2^(MixedGraphs::pa(x, v)-1))+1
  })

  ## set these to 1, add in vertex and set to -1
  out[pa_set] <- out[pa_set] + 1
  out[pa_set+2^(seq_len(n)-1)] <- out[pa_set+2^(seq_len(n)-1)] - 1

  out <- as.imset(out)

  out
}

##' @method standard_imset imset
##' @export
standard_imset.imset <- function(x) {
  n <- log2(length(x))
  out <- rev(subsetMatrix(n) %*% rev(1 - x))
  out <- as.imset(out)

  out
}

##' @method char_imset mixedgraph
##' @export
char_imset.mixedgraph <- function(x) {
  n <- length(x$vnames)
  out <- rep(0, 2^n)
  out[1] <- 1

  if (nedge(x, "undirected") > 0) {
    cl <- ADMGs2::cliques(x)
    comp <- unlist(lapply(cl, powerSet), recursive = FALSE)
    wh <- sapply(comp, function(v) sum(2^(v-1))+1)
    out[wh] <- 1
  }

  ht <- ADMGs2::headsTails(x, r = FALSE)
  if (length(ht$heads) > 0) {
    sets <- unlist(mapply(tailPowSet, ht$heads, ht$tails, SIMPLIFY = FALSE), recursive = FALSE)
    wh <- sapply(sets, function(v) sum(2^(v-1))+1)
    out[wh] = 1
  }

  out <- as.imset(out)

  out
}

##' Get characteristic imset from standard
##'
##' @param x a standard imset
##'
##'
##' @method char_imset imset
##' @export
char_imset.imset <- function(x) {
  n <- log2(length(x))

  out <- 1 - rev(abs(subsetMatrix(n)) %*% rev(x))
  out <- as.imset(out)

  out
}

##' @export
as.imset <- function(x) {
  n <- log2(length(x))
  if (n != round(n)) stop("Length not a power of 2")
  x <- as.integer(round(x))

  nms <- ""
  for (i in 1:n) {
    nms <- c(nms, paste(nms, i, sep=""))
  }
  nms[1] = "0"
  names(x) <- nms
  class(x) <- "imset"

  x
}

##' Print method for imsets
##' @export
print.imset <- function(x, only_nz = TRUE, ...) {
  n <- log2(length(x))
  cat("imset on ", n, " variables:\n", sep="")

  if (only_nz) {
    if (any(x != 0)) print.default(x[x != 0])
    else cat("(all entries zero)\n")
  }
  else print.default(x[seq_along(x)])
  invisible(x)
}

# imset2 <- function(graph) {
#   out <- rep(0,2^(graph$n))
#   out[1] = -1
#   out[length(out)] = 1
#
#   out2 <- out
#
#   tmp <- headsTails(graph)
#
#   wts <- 2^(seq(graph$n)-1)
#
#   for (i in seq_along(tmp$heads)) {
#     # sets <- set(tmp$heads[[i]], tmp$tails[[i]])
#     # locs <- wh_entries(sets)
#     # sgn2 <- setSign(tmp$heads[[i]], tmp$tails[[i]])
#     #
#     # out2[locs] <- out2[locs] + sgn2
#
#     sgn <- (-1)^(length(tmp$heads[[i]])-1)
#     small = sum(wts[tmp$tails[[i]]])
#     big = small + sum(wts[tmp$heads[[i]]])
#     out[small+1] = out[small+1] + sgn
#     out[big+1] = out[big+1] - sgn
#   }
#
#   # print(out-out2)
#   # print(out2)
#   out
# }

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
elemImset <- function(A, B, C=integer(0), n=max(c(A,B,C))) {
  out <- rep(0,2^n)

  wts <- 2^(seq_len(n)-1)
  wtC <- sum(wts[C])
  wtA <- sum(wts[A])
  wtB <- sum(wts[B])
  out[c(wtC, wtC+wtA, wtC+wtB, wtC+wtA+wtB)+1] = c(1,-1,-1,1)

  as.imset(out)
}

`+.imset` <- function(e1, e2) {
  if (length(e1) != length(e2)) stop("imsets defined over different subsets")
  NextMethod()
}

# gr <- graphCr("1->3<->2->4<->1", format = "ADMG")
# gr2 <- graphCr("1->3<-2->4<-1", format = "ADMG")
#
# all.equal(imset(gr), imset(gr2))
#
# imset(graphCr("1,2,3,4", format = "ADMG"))
