##' Compute Andrews' summation
##'
##' Compute the right-hand side of Andrews' entropy formula.
##'
##' @param graph a \code{mixedgraph} object of summary graph form
##' @param p a distribution or list or tables object of distributions
##'
##' @details This function uses the parameterizing sets to compute the entropy
##' with respect to each distribution in \code{p}.  If the distribution is in
##' the model then this will be equal to the entropy.
##'
##' @return A numeric vector with the same length as the number of distributions
##' supplied.
##'
##' @export
set_entropy <- function (graph, p) {

  n <- nv(graph)

  ## get parameterizing sets
  paramSets <- subsetRep(graph, sort = 3)
  ps <- powerSet(seq_len(n))

  ## determine which sets are in paramSets
  wh_set <- match(paramSets, ps)
  S <- subsetMatrix(n)

  ## now sum columns of subsetMatrix, and exclude values with a 0 coefficient
  coef <- colSums(S[wh_set,,drop=FALSE])
  sets <- ps[coef != 0]
  coef <- coef[coef != 0]

  set_entropy2(sets[-1], coef[-1], p)
}

## @param sets sets to use in entropy formula
## @param coef coefficients to apply
## @param p distribution, or a list or tables object of distributions
##' @importFrom contingency entropy
set_entropy2 <- function(sets, coef, p){
  if (!is.list(sets)) {
    return("input must be a list(set)")
  }
  if (length(sets) != length(coef)) stop("signs must have same length as sets")

  ## speed things up by removing zero coefficients
  if (any(coef == 0)) {
    sets <- sets[coef != 0]
    coef <- coef[coef != 0]
  }

  ## check input type
  if ("tables" %in% class(p)) out <- numeric(ntables(p))
  else if (is.list(p)) out <- numeric(length(p))
  else if (is.array(p) || is.numeric(p)) out <- numeric(1)
  else stop("Incorrect value at input")

  if (is.numeric(p)) {
    for (i in seq_along(sets)) {
      out <- out + coef[i]*entropy(p, sets[[i]])
    }
  }
  else if (is.list(p)) {
    out <- sapply(p, function(x) {
      for (i in seq_along(sets)) {
        out <- out + coef[i]*entropy(x, sets[[i]])
      }
    })
  }

  out
}
