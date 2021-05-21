##' @export
imset_conj1_v <- function(graph, v) {
  out <- elemImset(v, setdiff(graph$v, mb(graph, v)), setdiff(mb(graph, v), v))

  ## get heads and tails containing v
  ht <- headsTails(graph, r=FALSE)
  kp <- sapply(ht$heads, function(x) is.element(v,x))
  ht_v <- purrr::transpose(purrr::transpose(ht)[kp])

  for (h in seq_along(ht_v$heads)) {
    hd <- ht_v$heads[[h]]
    tl <- ht_v$tails[[h]]
    int <- ht_v$intrinsic[[h]]

    for (S in powerSet(setdiff(hd, v))) {
      gr2 <- graph[setdiff(c(hd, tl), S)]

      ## get factors for independence
      fct <- setdiff(mb(gr2, v, sort=2), v)
      oth <- setdiff(gr2$v, c(fct, v))

      out <- out + (-1)^(length(S)+1)*elemImset(v, oth, fct)
    }
  }

  out
}

##' @export
imset_conj1 <- function(graph) {
  topOrd <- topologicalOrder(graph)

  f <- function(k) {
    vs <- topOrd[seq_len(k)]
    imset_conj1_v(graph[vs], topOrd[k])
  }

  out <- zero_imset(nv(graph))
  for (i in seq_len(nv(graph))) out <- out + f(i)

  out
}

##' Recursive function to check if an imset if combinatorial
##'
##' @param imset imset to be checked
##'
##' @details This function recursively tries to remove
##'
##' @export
isCombinatorial <- function(imset) {
  return(isComb2(imset))
}

isComb2 <- function (imset) {
  if (all(imset == 0)) return(TRUE)
  wh <- max(which(imset != 0))
  if (imset[wh] < 0) return(FALSE)

  set <- posToSubset(wh, simplify = TRUE)
  if (length(set) < 2) return(FALSE)
  subsets <- combn(set, 2, simplify = FALSE)

  for (i in seq_along(subsets)) {
    tmp <- Recall(imset - elemImset(subsets[[i]][1], subsets[[i]][2],
                                    setdiff(set, subsets[[i]])))
    if (tmp) {
      print(c(subsets[[i]], setdiff(set, subsets[[i]])))
      return(TRUE)
    }
  }

  return(FALSE)
}
