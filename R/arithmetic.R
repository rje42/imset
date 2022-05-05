##' Overload arithmetic for imsets
##'
##' @param e1,e2 two imsets
##'
##' @details These functions are important, because they prevent accidentally
##' adding imsets with a different number of variables to one another, and the
##' recycling in R giving the wrong answer.
##'
##' @name imset_arith
NULL


##' @describeIn imset_arith Add imsets
##' @export
`+.imset` <- function(e1, e2) {

  if (length(e1) != length(e2)) {
    if (length(e1) > length(e2))
    {
      tmp <- e2
      e2 <- e1
      e1 <- tmp
    }
    out <- e2
    if (is.imset(e1)) out[seq_along(e1)] <- out[seq_along(e1)] + e1
    else out <- as.numeric(out) + e1

    return(as.imset(out))
  }

  NextMethod()
}

##' @describeIn imset_arith Subtract imsets
##' @export
`-.imset` <- function(e1, e2) {

  e2 <- (-1)*(e2)

  return(e1 + e2)
}
