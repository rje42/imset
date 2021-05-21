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

    for (S in powerSet(setdiff(hd, v))) {
      gr2 <- graph[setdiff(c(hd, tl), S)]

      ## get factors for independence
      fct <- setdiff(mb(gr2, v, sort=2), v)
      oth <- setdiff(gr2$v, c(fct, v))

      out <- out + (-1)^(length(hd)-length(S)+1)*elemImset(v, oth, fct)
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
