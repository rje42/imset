## This function attempts to split sets in a partition so that the set is rootable
##
## @param A the current partition
## @param Adj an adjacency matrix
## @param wh sequence representing the current splits
sort_partition <- function(A, Adj, wh=seq_along(A)) {

  ## get parameter and matrix to return
  lA <- length(A)
  M <- matrix(0, lA, lA)

  for (j1 in seq_len(lA)[-1]) {
    for (j2 in seq_len(j1-1)) {
      ## for each pair of entries in A, if they are all adjacent or not, then record this
      if (all(Adj[A[[j1]], A[[j2]]] == 1)) M[j1,j2] <- M[j2,j1] <- 1
      else if (all(Adj[A[[j1]], A[[j2]]] == 0)) M[j1,j2] <- M[j2,j1] <- 0
      else if (length(A[[j1]]) > 1 || length(A[[j2]]) > 1) {
        ## otherwise, figure out a split for the first set
        nmb1 <- rowSums(Adj[A[[j1]],A[[j2]],drop=FALSE])
        if (length(table(nmb1)) > 1 &&
            (all(nmb1 == 0 | nmb1 == length(j2)))) {
          A2 <- c(A[seq_len(j1-1)], list(A[[j1]][nmb1==0], A[[j1]][nmb1==length(j2)]), A[j1+seq_len(lA-j1)])
          wh <- wh[c(seq_len(j1),j1-1+seq_len(lA-j1+1))]
          return(Recall(A2, Adj, wh))
        }

        ## now try splitting the second set
        nmb2 <- colSums(Adj[A[[j1]],A[[j2]],drop=FALSE])
        if (length(table(nmb2)) > 1 &&
            (all(nmb2 == 0 | nmb2 == length(j1)))) {
          A2 <- c(A[seq_len(j2-1)], list(A[[j2]][nmb2==0], A[[j2]][nmb2==length(j1)]), A[j2+seq_len(lA-j2)])
          wh <- wh[c(seq_len(j2),j2-1+seq_len(lA-j2+1))]
          return(Recall(A2, Adj, wh))
        }

        ## note, no point trying to split on both sets
        ## if can't split on either individually
        ## so, if these both fail, then let the whole sequence fail
        return(NULL)
      }
      else {
        stop("We shouldn't get here")
      }
    }
  }

  attr(M, "seq") <- wh
  return(M)
}

##' @importFrom e1071 permutations

## check that variables in a partition are 'rooted' in the sense of Hu and Evans (2022)
is_rooted <- function(M) {

  ## if at most one set in the partition, there is nothing to show
  n <- nrow(M)
  if (n <= 1) return(TRUE)

  ## determine an ordering
  M2 <- 1*upper.tri(M)
  cs <- colSums(M == M2)
  cand_r <- which.max(cs) # get piece most connected to earlier pieces

  if (cs[cand_r] < n) {
    sq <- attr(M, "seq")
    if (any(duplicated(sq))) {
      tb <- table(sq)
      perms <- lapply(tb, function(x) e1071::permutations(x))
      k <- 0
      npms <- integer(length(tb))

      ## go through possible permutations check if any work
      for (i in seq_along(perms)) {
        d <- dim(perms[[i]])
        perms[[i]] <- (k + seq_len(ncol(perms[[i]])))[perms[[i]]]
        dim(perms[[i]]) <- d
        k <- k + d[2]
        npms[i] <- d[1]
      }
      combn <- combinations(npms)+1

      ## now go through combinations of permutations
      for (i in seq_len(nrow(combn))[-1]) {
        M2 <- M
        attr(M2, "seq") <- seq_along(ncol(M))
        prm <- integer(0)
        for (j in seq_along(perms)) {
          prm <- c(prm, perms[[j]][combn[i,j],])
        }
        ## check the partition for each one
        out <- Recall(M2[prm,prm])
        if (out) return(TRUE)
      }
    }
    return(FALSE)
  }
  if (!all(M[seq_len(cand_r-1), cand_r + seq_len(n - cand_r)] == 1)) {
    return(FALSE)
  }

  ## now split the model into the two pieces, and check each is also rooted
  Recall(M[seq_len(cand_r-1), seq_len(cand_r-1), drop=FALSE]) &&
    Recall(M[cand_r + seq_len(n - cand_r), cand_r + seq_len(n - cand_r), drop=FALSE])
}

checkDecomp <- function(graph, v) {
  n <- nv(graph)
  graph <- withAdjMatrix(graph)
  Adj <- collapse(graph$edges, dir=0L)

  ## if vertices not specified, check all of them
  if (missing(v)) v <- graph$v

  out <- rep(NA, max(v))

  ## go through vertices
  for (i in v) {
    fail <- FALSE
    nbs <- adj(graph, i)
    if (length(nbs) <= 1 || length(nbs) == n-1) {
      out[i] <- TRUE
      next
      # attr(out, "v") <- i
      # return(out)
    }

    ## get list of neighbours
    N <- vector(mode = "list", length = length(nbs))
    for (j in seq_along(nbs)) {
      N[[j]] <- setdiff(adj(graph, nbs[j]), c(nbs,i))
    }
    lens <- lengths(N)
    N <- N[order(lens)]
    nbs <- nbs[order(lens)]
    lens <- lens[order(lens)]
    for (j in seq_along(N)[-1]) {
      ## check neighbours follow sequential subsets
      if (!is.subset(N[[j-1]], N[[j]])) {
        fail <- TRUE
        out[i] <- FALSE
        break
      }
    }

    if (fail) next

    ## now check adjacencies within these neighbours
    N <- N[!duplicated(lens)]
    if (length(N) <= 1) {
      out[i] <- TRUE
      next
    }
    A <- vector("list", length=length(N))
    for (j in seq_along(A)) A[[j]] <- c(nbs[lens==unique(lens)[j]])

    ## implement the partition sorting function
    M <- sort_partition(A, Adj)
    if (is.null(M)) out[i] <- FALSE
    else out[i] <- is_rooted(M)

    # if (!fail) {
    #   out[i] <- TRUE
    #   # attr(out, "v") <- i
    #   # return(out)
    # }
    # else out[i] <- FALSE
  }

  return(out)
}

##' Implement algorithm from Hu and Evans (2022)
##'
##' @param graph bidirected graph of class \code{mixedgraph}
##'
##' @details Implements the algorithm described in Section 5 of Hu and Evans
##' (2022).
##'
##' @examples
##' data(bidi5)
##' bidi_decomposition(bidi5[[1]])
##' bidi_decomposition(bidi5[[1]])
##'
##' @export
bidi_decomposition <- function(graph) {
  ## check input is bidirected
  if (!is_BG(graph)) stop("'graph' must be a bidirected graph")

  if (nv(graph) <= 2) return(TRUE)
  ## if more than two vertices, check that the dual can be decomposed
  out1 <- checkDecomp(dual(graph))
  out2 <- logical(max(graph$v))
  for (v in which(out1)) {
    out2[v] <- Recall(graph[-v])
    if (out2[v]) return(TRUE)
  }
  return(any(out2, na.rm = TRUE))
}
