##' @export
stand_test <- function(graph, show_indeps=FALSE, signs=NULL) {
  v <- topologicalOrder(graph)

  ht <- headsTails(graph, r=FALSE)
  max_v <- sapply(ht$heads, function(x) max(match(x, v)))
  ord <- order(max_v)
  ht <- list(heads=ht$heads[ord], tails=ht$tails[ord], intrinsic=ht$intrinsic[ord])
  class(ht) <- "htList"
  max_v <- v[max_v[ord]]

  out <- zero_imset(length(graph$vnames))

  for (i in seq_along(v)[-1]) {
    pv <- v[seq_len(i-1)]
    mbl <- mb(graph, v[i], c(pv,v[i]))
    RHS <- setdiff(pv, mbl)
    if (show_indeps && length(RHS) > 0) print(as.ci(v[i], RHS, setdiff(mbl, v[i])))
    out <- out + elem_imset(v[i], RHS, setdiff(mbl, v[i]))

    rel_heads <- which(max_v == v[i])

    ## now add in head independences
    for (j in seq_along(rel_heads)) {
      hd <- ht$heads[rel_heads][[j]]
      tl <- ht$tails[rel_heads][[j]]
      if (isTRUE(all.equal(hd, v[i]))) next

      anH <- anc(graph, hd)
      hd_i <- setdiff(hd, v[i])

      if (is.null(signs)) sgn <- -kronPower(c(1,-1), length(hd)-1)[-1]
      else {
        sgn <- signs
        sgn <- rep_len(sgn, 2^(length(hd)-1)-1)
      }
      ps <- powerSet(hd_i)[-1]
      mbs <- lapply(ps, function(x) mb(graph, v[i], setdiff(anH, x)))
      sol <- mapply(function(x,y,z) z*elem_imset(v[i], setdiff(c(hd_i,tl),c(x,y)), setdiff(x,v[i])), mbs, ps, sgn, SIMPLIFY = FALSE)
      if (length(sgn) != length(ps)) stop("Problem with sgn")
      if (show_indeps) {
        cat("head: ", paste(hd, collapse=", "), "\n")
        for (k in seq_along(sgn)) {
          RHS <- setdiff(c(hd_i,tl),c(mbs[[k]],ps[[k]]))
          cond <- setdiff(mbs[[k]], v[i])
          if (length(RHS) > 0) cat(sgn[k], "  ", v[i], "_||_", RHS, ifelse(length(cond) > 0, "|", ""), cond, "\n")
        }
      }

      out <- out + Reduce(`+`, sol)
    }

  }
  out
}

##' @export
si_district <- function (graph) {
  dist <- districts(graph)
  dist <- dist[[which.max(lengths(dist))]]
  not_v <- setdiff(graph$v, dist)
  for (i in not_v) graph <- mutilate(graph, i, etype=c("directed", "bidirected"), dir=c(-1L,0L))

  if (length(not_v) > 1) {
    ne <- combn(not_v, 2)
    class(ne) <- "edgeMatrix"
    graph <- addEdges(graph, list(dir=ne))
  }

  standard_imset(graph)
}
