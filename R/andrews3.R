##' Convert between ci object and text representation
##'
##' @param cis list of ci objects
##' @param symm logical: should symmetric independences be listed
##'
##' @export
ci2andrews <- function (cis, symm=FALSE) {
  if (length(cis) == 0) {
    if (symm) return(matrix(character(0), nrow=2))
    else return(character(0))
  }
  if (!is.list(cis[[1]])) cis <- list(cis)

  out <- sapply(cis,
                function (ci) if (length(ci[[3]]) > 0) paste(paste0(ci[[1]], collapse=""), ":",
                                                             paste0(ci[[2]], collapse=""), "|",
                                                             paste0(ci[[3]], collapse=""), sep="")
                else paste(paste0(ci[[1]], collapse=""), ":",
                           paste0(ci[[2]], collapse=""), sep=""))
  if (symm) {
    ## add in symmetric independences if required
    out2 <- sapply(cis,
                   function (ci) if (length(ci[[3]]) > 0) paste(paste0(ci[[2]], collapse=""), ":",
                                                                paste0(ci[[1]], collapse=""), "|",
                                                                paste0(ci[[3]], collapse=""), sep="")
                   else paste(paste0(ci[[2]], collapse=""), ":",
                              paste0(ci[[1]], collapse=""), sep=""))
    out <- c(rbind(out, out2))
  }
  return(out)
}

##' @describeIn ci2andrews convert text string to ci objects
##' @param x string of Andrew's notation
##' @export
andrews2ci <- function (x) {
  ## start with some string manipulation
  out <- strsplit(x, split=",")[[1]]
  pval <- regexpr("[|]", out)
  colval <- regexpr("[:]", out)
  end <- nchar(out)
  cond <- (pval > 0)

  ## now deduce triples
  A <- B <- C <- vector(mode="list", length=length(out))
  A <- lapply(strsplit(sapply(strsplit(out, split=":"), function(x) x[[1]]), ""), as.numeric)
  B[cond] <- mapply(function(x,y,z) substr(x,y,z), out[cond], y=colval[cond]+1, z=pval[cond]-1)
  B[!cond] <- mapply(function(x,y,z) substr(x,y,z), out[!cond], y=colval[!cond]+1, z=end[!cond])
  B <- lapply(unlist(B), function(x) strsplit(x, "")[[1]])
  C[cond] <- mapply(function(x,y,z) substr(x,y,z), out[cond], y=pval[cond]+1, z=end[cond])
  C[cond] <- lapply(unlist(C[cond]), function(x) strsplit(x, "")[[1]])
  C[!cond] <- list(integer(0))

  ## then return the solution
  mapply(as.ci, A, B, C, SIMPLIFY = FALSE)
}

# NIE <- function (graph, topOrd) {
#   if (missing(topOrd)) topOrd <- topologicalOrder(graph)
#
#   incL <- excL <- list()
#
#   out <- pairs(graph, topOrd=topOrd)
#   n <- length(topOrd)
#
#   ## now construct imset
#   for (i in rev(seq_len(n))) {
#     if (length(out[[i]]$N) == 0) next
#
#     ps <- powerSet(seq_len(i))
#     NJ <- MK <- list(topOrd[seq_len(i)])
#     for (J in seq_along(ps)[-1]) {
#       NJ[[J]] <- Reduce(intersect, out[[i]]$N[ps[[J]]])
#       MK[[J]] <- Reduce(intersect, out[[i]]$M[ps[[J]]])
#     }
#
#     subM <- subsetMatrix(i)
#
#     if (i < length(topOrd)) ext <- elem_imset(topOrd[seq_len(n-i)+i],
#                                               topOrd[seq_len(i)])
#     else ext <- zero_imset(n)
#
#     for (j in seq_along(NJ)) for (k in which(subM[j,] != 0)) {
#       B <- setdiff(NJ[[j]], MK[[k]])
#       C <- setdiff(intersect(NJ[[j]], MK[[k]]), topOrd[i])
#       # cat(subM[j,k], ": ")
#       # print(list(topOrd[i], B, C))
#       # print(1 - char_imset(elem_imset(topOrd[i], B, C)))
#       sig <- list(c(2^(topOrd[i]-1), sum(2^(B-1)), sum(2^(C-1))))
#       if (subM[j,k] > 0) {
#         incL <- c(incL, sig)
#       }
#       else if (subM[j,k] < 0) {
#         excL <- c(excL, sig)
#       }
#       else stop("Shouldn't get to here")
#     }
#   }
#
#   inc <- exc <- zero_imset(nv(graph))
#   for (i in seq_along(incL)) {
#     sts <- int2set(incL[[i]])
#     inc <- inc + 1 - char_imset(elem_imset(sts[[1]], sts[[2]], sts[[3]], n=n))
#   }
#   for (i in seq_along(excL)) {
#     sts <- int2set(excL[[i]])
#     exc <- exc + 1 - char_imset(elem_imset(sts[[1]], sts[[2]], sts[[3]], n=n))
#   }
#
#   return(list(inc=inc, exc=exc))
# }

##' Implement Andrews' Algorithm 6
##'
##' @param graph an ADMG of class \code{mixedgraph}
##' @param topOrd an optional topological ordering
##' @param alg3 logical: should the less sophisticated Algorithm 3 be used?
##'
##' @details Implements Algorithm 6 from Andrews (2022).  This involves, for
##' each vertex v, first collecting pairs of maximal disconnected sets and the
##' corresponding Markov blankets of v, and then taking all intersections of
##' these sets.  Each of these can be associated with an independence, which is
##' either added to the inclusion or the exclusion set, depending upon its
##' parity.  The difference between these two imsets is the characteristic imset
##' that can be obtained directly using \code{char_imset()}.
##'
##' Algorithm 3 (which is run if \code{alg3} is \code{TRUE}) is essentially
##' the same as Algorithm 6, but it does not check for independences that have
##' already been added so ends up with much higher degree results for both the
##' inclusion and exclusion imsets.  The difference between the two is the same
##' for both, however.
##'
##' @export
NIE <- function (graph, topOrd, alg3 = FALSE) {
  if (missing(topOrd)) topOrd <- topologicalOrder(graph)

  incL <- excL <- list()

  out <- pairs(graph, topOrd=topOrd)
  n <- length(topOrd)

  ## now construct imset
  for (i in rev(seq_len(n))) {
    if (length(out[[i]]$N) == 0) next

    ps <- powerSet(seq_len(i))
    NJ <- MK <- list(topOrd[seq_len(i)])
    for (J in seq_along(ps)[-1]) {
      NJ[[J]] <- Reduce(intersect, out[[i]]$N[ps[[J]]])
      MK[[J]] <- Reduce(intersect, out[[i]]$M[ps[[J]]])
    }

    subM <- subsetMatrix(i)

    for (j in seq_along(NJ)) for (k in which(subM[j,] != 0)) {
      B <- setdiff(NJ[[j]], MK[[k]])
      C <- setdiff(intersect(NJ[[j]], MK[[k]]), topOrd[i])
      # cat(subM[j,k], ": ")
      # print(list(topOrd[i], B, C))
      # print(1 - char_imset(elem_imset(topOrd[i], B, C)))
      sig <- list(c(2^(topOrd[i]-1), sum(2^(B-1)), sum(2^(C-1))))
      if (subM[j,k] > 0) {
        if (!alg3 && match(sig, excL, nomatch = 0L) > 0) excL <- excL[-match(sig, excL, nomatch = 0L)]
        else incL <- c(incL, sig)
      }
      else if (subM[j,k] < 0) {
        if (!alg3 && match(sig, incL, nomatch = 0L) > 0) incL <- incL[-match(sig, incL, nomatch = 0L)]
        else excL <- c(excL, sig)
      }
      else stop("Shouldn't get to here")
    }
  }

  inc <- exc <- zero_imset(nv(graph))
  for (i in seq_along(incL)) {
    sts <- int2set(incL[[i]])
    inc <- inc + 1 - char_imset(elem_imset(sts[[1]], sts[[2]], sts[[3]], n=n))
  }
  for (i in seq_along(excL)) {
    sts <- int2set(excL[[i]])
    exc <- exc + 1 - char_imset(elem_imset(sts[[1]], sts[[2]], sts[[3]], n=n))
  }

  return(list(inc=inc, exc=exc))
}
