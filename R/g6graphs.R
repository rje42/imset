collect_graph <- function(i, file="6c", dir="graphs", os="mac") {
  cmd <- paste0("sed '", i, "q;d' ", dir, "/graph", file, ".g6 | ~/showg_", os, "64")

  tmp <- system(cmd, intern=TRUE)
  if (length(tmp) == 0) stop("Beyond limits of file")
  tmp <- tmp[-(1:2)]
  tmp <- strsplit(tmp, " : ")
  tmp <- sapply(tmp, function(x) x[2])
  edges <- lapply(tmp, function(x) as.numeric(strsplit(substr(x, 1, nchar(x)-1), " ")[[1]])+1)

  ## now we have adjList, so create graph
  class(edges) <- "adjList"
  edg <- list(undirected=edges)
  class(edg) <- "edgeList"

  return(mixedgraph(n=length(edges), edges=edg))
}
