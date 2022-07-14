##' Turn conditional independence
ci2entropy <- function(...) {
  UseMethod("ci2entropy", ...)
}

ci2entropy.ci <- function(ci, n) {
  if (missing(n)) n <- max(unlist(ci))
  ci2entropy.default(ci[[1]], ci[[2]], ci[[3]], n=n)
}

ci2entropy.default <- function(A, B, C, n) {
  if (missing(n)) n <- max(c(A,B,C))

  if (missing(C)) C = integer(0)
  else {
    if (length(intersect(A,C) > 0)) A <- setdiff(A,C)
    if (length(intersect(B,C) > 0)) B <- setdiff(B,C)
  }
  if (length(intersect(A,B)) > 0) stop("Sets A and B should not intersect")

  (-1)*elem_imset.default(A=A, B=B, C=C, n=n)
}

# ##' Shannon Cone
# shannon_cone <- function(n) {
#   -gen_constraints(n)
# }
