##' Recursive function to check if an imset if combinatorial
##'
##' @param imset imset to be checked
##' @param use_char logical: should characteristic imset method be used?
##' @param trace logical: show detail
##'
##' @details This function recursively tries to remove conditional independences
##' by greedy searching.  It looks at the largest set still present, and then
##' tries all pairs still in the characteristic imset that are contained in
##' this set as variables to be made conditionally independent.  The associated
##' elementary characteristic imset is then subtracted, until we either reach
##' the zero imset (in which case the imset is combinatorial), or we see a
##' negative value.  In the latter case, we reject that sequence and try
##' another branch.  Once all branches are exhausted the algorithm fails.
##'
##' @seealso \code{\link{defines_mod}}
##'
##' @export
is_combinatorial <- function(imset, use_char=TRUE, trace=FALSE) {
  if (use_char) {
    imset2 <- 1 - char_imset.imset(imset)
    out <- isComb3(imset2, trace=trace)
  }
  else out <- isComb2(imset, trace=trace)
  class(out) <- "isComb"
  return(out)
}


##' Deprecated recursive function to check if an imset if combinatorial
##'
##' @param imset imset to be checked
##' @param use_char logical: should characteristic imset method be used?
##' @param trace logical: show detail
##'
##' @details This function recursively tries to remove conditional independences
##' by greedy searching.  It looks at the largest set still present, and then
##' tries all pairs still in the characteristic imset that are contained in
##' this set as variables to be made conditionally independent.  The associated
##' elementary characteristic imset is then subtracted, until we either reach
##' the zero imset (in which case the imset is combinatorial), or we see a
##' negative value.  In the latter case, we reject that sequence and try
##' another branch.  Once all branches are exhausted the algorithm fails.
##'
##' It has been replaced by \code{is_combinatorial}
##'
##' @seealso \code{\link{defines_mod}}, \code{\link{is_combinatorial}}
##'
##' @name isCombinatorial-deprecated
##' @export
isCombinatorial <- function(imset, use_char=TRUE, trace=FALSE) {
  .Deprecated(is_combinatorial, msg="This function is deprecated, use is_combinatorial instead")
  cl <- match.call()
  args <- as.list(cl[-1])
  do.call(is_combinatorial, args)
}

isComb2 <- function (imset, trace=FALSE) {
  if (trace > 0) cat(paste(paste(rep("-", trace-1), collapse=""), "check values... "))
  if (all(imset == 0)) {
    if (trace > 0) cat("\n")
    return(TRUE)
  }
  wh <- max(which(imset != 0))
  if (imset[wh] < 0) {
    if (trace > 0) cat("\n")
    return(FALSE)
  }

  set <- posToSubset(wh, simplify = TRUE)
  if (length(set) < 2) {
    if (trace > 0) cat("\n")
    return(FALSE)
  }
  subsets <- combn(set, 2, simplify = FALSE)

  for (i in seq_along(subsets)) {
    if (trace > 0) {
      cat(paste("trying ", subsets[[i]][1],",", subsets[[i]][2],
                " | ", paste(setdiff(set, subsets[[i]]),collapse=""), "\n", sep=""))
      # cat(paste("trying ", subsets[[i]][1],",", subsets[[i]][2],
      #           paste(setdiff(set, subsets[[i]]),collapse=""), "\n"))
      trace = trace + 1
    }
    tmp <- Recall(imset - elem_imset(subsets[[i]][1], subsets[[i]][2],
                                    setdiff(set, subsets[[i]])), trace=trace)
    if (tmp) {
      if (trace > 0) print(c(subsets[[i]], setdiff(set, subsets[[i]])))
      out <- TRUE
      attr(out, "elem") <- c(attr(tmp, "elem"), list(elem_imset(subsets[[i]][1], subsets[[i]][2],
                                                               setdiff(set, subsets[[i]]))))
      return(out)
    }
  }
  if (trace > 0) cat("path failed\n")
  return(FALSE)
}

isComb3 <- function (imset, trace=FALSE) {
  if (trace > 0) cat(paste(paste(rep("-", trace-1), collapse=""), "check values... "))
  if (all(imset == 0)) {
    if (trace > 0) cat("\n")
    return(TRUE)
  }
  if (any(imset < 0)) {
    if (trace > 0) cat("\n")
    return(FALSE)
  }

  wh <- max(which(imset > 0))

  set <- posToSubset(wh, simplify = TRUE)
  if (length(set) < 2) {
    if (trace > 0) cat("\n")
    return(FALSE)
  }

  ## get possible pairs to check
  v <- as.integer(intToBits(sum(2^(set-1))))
  comp <- sapply(which(imset > 0), function(x) as.integer(intToBits(x-1)))
  if (is.null(dim(comp))) dim(comp) <- c(32, length(comp)/32)

  colsKp <- colSums(comp) == 2 & apply(v - comp, 2, function(x) all(x >= 0))
  comp <- comp[,colsKp,drop=FALSE]
  if (ncol(comp) == 0) {
    if (trace > 0) cat("\n")
    return(FALSE)
  }

  for (i in seq_len(ncol(comp))) {
    x <- which(comp[,i] > 0)
    y <- x[2]
    x <- x[1]

    if (trace > 0) {
      cat(paste(paste(rep("-", trace-1), collapse=""), " trying ", x,",", y,
                " | ", paste(setdiff(set, c(x,y)),collapse=""), "\n", sep=""))
      trace = trace + 1
    }
    tmp <- Recall(imset - 1 + char_imset(elem_imset(x, y, setdiff(set, c(x,y)),
                                                   n = log2(length(imset)))),
                  trace=trace)
    if (tmp) {
      if (trace > 0) print(c(x, y, setdiff(set, c(x,y))))
      out <- TRUE
      attr(out, "elem") <- c(attr(tmp, "elem"), list(elem_imset(x, y,
                                                               setdiff(set, c(x,y)))))
      return(out)
    }
  }
  if (trace > 0) cat("path failed\n")
  return(FALSE)
}

##' @exportS3Method print isComb
print.isComb <- function(x, ...) {
  if (isFALSE(x)) print(FALSE)
  else {
    print(TRUE)
    CIs <- lapply(attr(x, "elem"), sei2ci)
    if (length(CIs) > 0) {
      cat("This is a combinatorial imset, and is made up of the following independences: \n")
      for (i in seq_along(CIs)) {
        print(CIs[[i]])
      }
    }
    else cat("This is an empty and (therefore) combinatorial imset\n")
  }

  invisible(x)
}
