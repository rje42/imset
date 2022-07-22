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
##'
##' @details Implements Algorithm 6 from Andrews (2022).
##'
##' @export
NIE <- function (graph, topOrd) {
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
        if (match(sig, excL, nomatch = 0L) > 0) excL <- excL[-match(sig, excL, nomatch = 0L)]
        else incL <- c(incL, sig)
      }
      else if (subM[j,k] < 0) {
        if (match(sig, incL, nomatch = 0L) > 0) incL <- incL[-match(sig, incL, nomatch = 0L)]
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
