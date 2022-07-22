pairs <- function (graph, topOrd, subset_rep) {
  if (missing(topOrd)) topOrd <- topologicalOrder(graph)
  if (missing(subset_rep)) subset_rep <- subsetRep(graph, sort=3)

  ord <- sapply(subset_rep, function (x) max(match(x, topOrd)))
  subs <- vector(length=nv(graph), mode="list")
  for (i in seq_along(subs)) {
    subs[[i]] <- subset_rep[ord == topOrd[i]]
  }

  out <- vector(length=nv(graph), mode="list")
  names(out) <- graph$vnames[topOrd]
  for (i in seq_along(topOrd)) {
    out[[i]] <- pairs_v(graph[topOrd[seq_len(i)]], topOrd[i], subs[[i]])
  }

  out
}

pairs_v <- function (graph, v, subset_rep) {

  # if (missing(topOrd)) topOrd <- topologicalOrder(graph)

  ## get position in topological ordering
  # vpos <- match(v, topOrd)
  #
  # if (!missing(subset_rep)) {
  #   subset_rep <- subset_rep[sapply(subset_rep, function(x) !all(x %in% topOrd[seq_len(vpos)]))]
  # }

  ## determine maximal vertex
  # graph <- graph[topOrd[seq_len(vpos)], drop=TRUE]

  # if (missing(subset_rep)) subset_rep <- subsetRep(graph, sort = 3)

  sub_mis <- powerSetCond(v, graph$v[-length(graph$v)], sort = 3)
  sub_mis <- sub_mis[!(sub_mis %in% subset_rep)]

  lm <- lengths(subset_rep)
  ln <- lengths(sub_mis)

  i <- 0
  N <- M <- list()

  while (length(sub_mis) > 0) {
    i <- i+1
    N[[i]] <- sub_mis[[which.max(ln)]]
    subM <- sapply(subset_rep, is.subset, y=N[[i]])
    # Ms <- sapply(subset_rep, is.subset, y=N[[i]])
    M[[i]] <- subset_rep[[which.max(lm*subM)]]
    subM <- sapply(sub_mis, is.subset, y=N[[i]]) & !sapply(sub_mis, is.subset, y=M[[i]])
    sub_mis <- sub_mis[!subM]
    ln <- ln[!subM]
  }

  list(N=N, M=M)
}

