##' Construct sparse imset
##'
##' @param x entries
##' @param pos indices or list of subsets
##' @param n number of variables
##' @param vnames names for variables
##' @param names logical: should variable names be included?
##'
##'
make_imset <- function (x, pos, n, vnames, names=TRUE) {
  if (length(x) != length(pos)) stop("Entries and indices must have same length")

  if (is.list(pos)) {
    posi <- subsetToPos(pos)
  }
  else if (is.numeric(pos)) posi <- pos
  if (missing(n) && missing(vnames)) {
    n <- ceiling(log2(max(posi)))
  }
  else if (missing(n)) n <- length(vnames)

  if (names && missing(vnames)) vnames <- paste("x", seq_len(n), sep="")

  if (any(posi > 2^n)) stop("n provided not compatible with some entries")

  if (names) {
    if (is.list(pos)) names(x) <- entry_names.list(pos, vnames)
    else names(x) <- entry_names.numeric(posi, vnames)
  }

  ## sort entries
  ord <- order(posi)
  x <- x[ord]
  posi <- posi[ord]

  class(x) <- "sparse_imset"
  attr(x, "n") <- n
  attr(x, "idx") <- posi

  return(x)
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

##' @export
`+.sparse_imset` <- function (e1, e2) {

  p1 <- attr(e1, "idx"); p2 <- attr(e2, "idx")
  n1 <- attr(e1, "n"); n2 <- attr(e2, "n")
  l1 <- length(e1); l2 <- length(e2)
  n <- max(n1, n2)
  # m <- min(n1, n2)
  # if (!isTRUE(all.equal(e1$vnames[seq_len(m)], e2$vnames[seq_len(m)]))) stop("Variable names do not match")

  # ##
  # if (n1 >= n2) vnms <- e1$vnames
  # else vnms <- e2$vnames

  x <- posi <- integer(0)
  j <- 1
  j1 <- j2 <- 1

  ## go through the entries adding them appropriately
  while (j1 <= l1 && j2 <= l2) {
    if (p1[j1] > p2[j2]) {
      x[j] <- e2[j2]
      posi[j] <- p2[j2]
      j2 <- j2 + 1
      j <- j + 1
    }
    else if (p1[j1] < p2[j2]) {
      x[j] <- e1[j1]
      posi[j] <- p1[j1]
      j1 <- j1 + 1
      j <- j + 1
    }
    else if (e1[j1] + e2[j2] != 0) {
      x[j] <- e1[j1] + e2[j2]
      posi[j] <- p2[j2]
      j1 <- j1 + 1
      j2 <- j2 + 1
      j <- j + 1
    }
    else {
      j1 <- j1 + 1
      j2 <- j2 + 1
    }
  }

  ## now take tail of longer imset, if any
  if (j1 > l1 && j2 <= l2) {
    x <- c(x, e2[seq(j2,l2)])
    posi <- c(posi, p2[seq(j2,l2)])
  }
  else if (j2 > l2 && j1 <= l1) {
    x <- c(x, e1[seq(j1,l1)])
    posi <- c(posi, p1[seq(j1,l1)])
  }

  ## get resulting imset
  out <- make_imset(x, posi, n=n, names=!is.null(names(e1)))
  return(out)
}

elem_imset_sparse <- function(A, B, C = integer(0), n = max(c(A, B, C))) {
  return(make_imset(x = c(1,-1,-1,1), pos = list(c(A,B,C),c(B,C),c(A,C),C), n=n, names=FALSE))
}
