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
  out[subsetToPos(list(A))] <- 1

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
standard_imset.mixedgraph <- function(x, slow=FALSE) {
  if (!is_ADMG(x)) stop("This function only currently works for ADMGs")

  n <- length(x$vnames)
  out <- rep(0, 2^n)
  max_ent <- sum(2^(x$v-1))+1
  out[1] = -1; out[max_ent] = 1

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
  if (is_ADMG(x)) {
    if (!slow) ht <- headsTails3(x, r=FALSE, sort=3)
    else ht <- headsTails(x, r=FALSE, sort=3)

    for (h in seq_along(ht$heads)) {
      idx <- sapply(powerSet(ht$heads[[h]]), function(x) sum(2^(x-1))) + 1
      idx <- idx + sum(2^(ht$tails[[h]]-1))
      sgn <- kronPower(c(-1,1), length(ht$heads[[h]]))
      out[idx] <- out[idx] - sgn
    }
  }
  else if (is_UG(x)) {
    clq <- cliques(x)
    stop("Undirected graphs not yet implemented")
  }
  else stop("Graph must be an ADMG")

  out <- as.imset(out)

  out
}

##' @method standard_imset imset
##' @export
standard_imset.imset <- function(x) {
  n <- log2(length(x))
  out <- rev(invMobius(rev(1 - x)))
  # out <- rev(subsetMatrix(n) %*% rev(1 - x))
  out <- as.imset(out)

  out
}

##' @method char_imset mixedgraph
##' @export
char_imset.mixedgraph <- function(x) {
  if (!is_ADMG(x)) stop("This function only currently works for ADMGs")
  n <- length(x$vnames)
  out <- rep(0, 2^n)
  out[1] <- 1

  if (nedge(x, "undirected") > 0) {
    cl <- ADMGs2::cliques(x)
    comp <- unlist(lapply(cl, powerSet), recursive = FALSE)
    wh <- sapply(comp, function(v) sum(2^(v-1))+1)
    out[wh] <- 1
  }

  ht <- ADMGs2::headsTails3(x, r = FALSE)
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

  out <- 1 - rev(fastMobius(rev(x)))
  # out <- 1 - rev(abs(subsetMatrix(n)) %*% rev(x))
  out <- as.imset(out)

  out
}

##' Get maximal positive/negative sets from imset
##'
##' Intended for a structural imset, obtain the maximal positive and negative
##' sets in an imset.
##'
##' @param x imset
##' @name envelopes
NULL

##' @describeIn envelopes Get maximal positive sets from imset
##' @export
upperEnvelope <- function(x) {
  n <- log2(length(x))
  Sm <- subsetMatrix(n)
  wh <- which(x != 0)

  out <- list()

  while(any(x > 0)) {
    if (x[last(wh)] < 0) stop("Not a stuctural imset")
    out <- c(out, posToSubset(last(wh)))

    ## now set subsets of this one to zero
    x[Sm[last(wh),] != 0] <- 0
    wh <- which(x != 0)
  }
  out <- unname(out)
  out
}

## Get maximal negative sets from imset
##' @describeIn envelopes Get maximal negative sets from imset
##' @export
lowerEnvelope <- function(x) {
  n <- log2(length(x))
  Sm <- subsetMatrix(n)
  wh <- which(x < 0)

  out <- list()

  while(any(x < 0)) {
    out <- c(out, posToSubset(last(wh)))

    ## now set subsets of this one to zero
    x[Sm[last(wh),] != 0] <- 0
    wh <- which(x < 0)
  }
  out <- unname(out)
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
#     # locs <- setsToPos(sets)
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

