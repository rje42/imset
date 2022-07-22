##' Get the imset representing constrained sets from a conditional independence
##'
##' @param ci a conditional independence as an object of class \code{ci}
##' @param n total number of variables
##'
##' @export
cons_imset <- function(ci, n) {
  if (missing(n)) {

    n <- tryCatch(max(unlist(ci)), warning=0L)
  }

  out <- zero_imset(n)
  pos <- 2^(seq_len(n)-1)
  posList <- vector(mode="list", length=3)

  ## get positions in imset
  for (i in 1:3) posList[[i]] <- c(combinations(rep(2,length(ci[[i]]))) %*% pos[ci[[i]]])
  posList[[1]] <- posList[[1]][-1]
  posList[[2]] <- posList[[2]][-1]

  out[kronecker(kronecker(posList[[1]], posList[[2]], "+"),
                posList[[3]], "+") + 1] <- 1

  out
}

cons_sets <- function(imset) {
  if (any(imset < 0)) {
    warning("Negative entries, so converting to 1 - characteristic imset")
    imset <- 1 - char_imset(imset)
    if (any(imset < 0)) stop("Not a valid imset to find subsets")
  }

  wh <- int2set(which(imset > 0) - 1)
  return(wh)
}
