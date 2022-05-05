##' Test and enforce that vector is an imset
##'
##' @param x vector test or make imset from
##'
##' @details \code{as.imset} forces a vector to become an imset;
##' the length of the vector \code{x} must be a power of 2 length.
##' \code{is.imset} just checks that \code{imset} is in
##' \code{class(x)}.
##'
##' @name imset_control
NULL

##' @describeIn imset_control force object to be an imset
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

##' @describeIn imset_control test if object is an imset
##' @export
is.imset <- function(x) {
  return("imset" %in% class(x))
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
