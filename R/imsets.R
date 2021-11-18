# library(MixedGraphs)
# library(ADMGs)


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

##' Zero imset
##'
##' @param n number of variables
##'
##' @export
zero_imset <- function(n) {
  as.imset(rep(0, 2^n))
}

##' @describeIn zero_imset identify a single set
##' @param A set to identify
##' @export
identifier_imset <- function(A, n) {
  if (missing(n)) n <- max(A)

  out <- as.imset(rep(0, 2^n))
  out[wh_entries(list(A))] <- 1

  out
}

##' Standard imset
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
##' @importFrom rje kronPower
##' @export
standard_imset.mixedgraph <- function(x) {
  if (!is.ADMG(x)) stop("This function only currently works for ADMGs")

  n <- length(x$vnames)
  out <- rep(0, 2^n)
  out[1] = -1; out[2^n] = 1

  # if (is.DAG(x)) {
  #   ## get entries associated with parent sets
  #   pa_set <- sapply(seq_len(n), function(v) {
  #     sum(2^(MixedGraphs::pa(x, v)-1))+1
  #   })
  #
  #   pas <- table(pa_set)
  #   for (i in seq_along(pas)) {
  #     out[as.numeric(names(pas)[i])] <- out[as.numeric(names(pas)[i])] + pas[i]
  #   }
  #
  #   ## set these to 1, add in vertex and set to -1
  #   out[pa_set] <- out[pa_set] + 1
  #   out[pa_set+2^(seq_len(n)-1)] <- out[pa_set+2^(seq_len(n)-1)] - 1
  # }
  {
    ht <- headsTails(x, r=FALSE)

    for (h in seq_along(ht$heads)) {
      idx <- sapply(powerSet(ht$heads[[h]]), function(x) sum(2^(x-1))) + 1
      idx <- idx + sum(2^(ht$tails[[h]]-1))
      sgn <- kronPower(c(-1,1), length(ht$heads[[h]]))
      out[idx] <- out[idx] - sgn
    }
  }

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
##' @export
elemImset <- function(A, B, C=integer(0), n=max(c(A,B,C))) {

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

##' @describeIn elemImset Construct conditional independence from elementary imset
##' @export
elemToIndep <- function(imset) {
  if (sum(imset != 0) != 4) stop("Not an elementary imset")

  ## isolate sets
  n <- log2(length(imset))
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
