##' Construct sparse imset
##'
##' @param x entries
##' @param pos indices or list of subsets
##' @param n number of variables
##' @param vnames names for variables
##'
##'
make_imset <- function (x, pos, n, vnames) {
  if (length(x) != length(pos)) stop("Entries and indices must have same length")

  if (is.list(pos)) {
    posi <- subsetToPos(pos)
  }
  else if (is.numeric(pos)) posi <- pos
  if (missing(n) && missing(vnames)) {
    n <- ceiling(log2(max(posi)))
  }
  else if (missing(n)) n <- length(vnames)

  if (missing(vnames)) vnames <- paste("x", seq_len(n), sep="")

  if (any(posi > 2^n)) stop("n provided not compatible with some entries")

  if (is.list(pos)) names(x) <- entry_names.list(pos, vnames)
  else names(x) <- entry_names.numeric(posi, vnames)

  class(x) <- "sparse_imset"
  attr(x, "n") <- n
  attr(x, "idx") <- posi

  x
}

##' @exportS3Method print sparse_imset
print.sparse_imset <- function (x) {
  cat(sprintf("Sparse imset on %d variables:\n", attr(x, "n")))
  if (sum(x != 0) == 0) cat("(all entries zero)")
  else print.default(x[x != 0])
}

identifier_sparse_imset <- function(A, n, vnames) {
  make_imset(x=1, pos=subsetToPos(list(A)), n=n, vnames)
}

# `+.sparse_imset` <- function (e1, e2) {
#   pos1 <- attr(e1, "idx")
#   pos2 <- attr(e2, "idx")
#
#   out_i <- c(pos1, pos2)
#   ord <- order(out_i)
#
#   if ()
# }
