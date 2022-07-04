gen_constraints <- function (n, sparse=FALSE, reduce=NULL) {
  if (sparse) require(Matrix)
  if (n < 0) stop("n must be non-negative")
  if (n <= 1) return(matrix(NA, ncol=0, nrow=n+1))

  ## get pairs
  prs <- combn(n, 2) - 1
  idx <- 0
  if (!is.null(reduce)) reduce <- lapply(reduce, sort.int)
  else reduce <- list()

  if (sparse) {
    M <- new("dgTMatrix",
             i = integer(0),
             j = integer(0),
             x = numeric(0),
             Dim = as.integer(c(2^n, choose(n,2)*2^(n-2))))
  }
  else M <- matrix(0, 2^n, choose(n,2)*2^(n-2))

  colnames(M) <- character(ncol(M))

  for (i in seq_len(choose(n,2))) {
    if (setmatch(list(prs[,i]+1), reduce, nomatch = 0L) > 0) next
    pr <- prs[,i]
    vals <- c(0, 2^pr[1], 2^pr[2], 2^pr[1] + 2^pr[2])

    kp <- rep(c(TRUE,FALSE), each=2^pr[1])
    kp <- kp*rep(c(TRUE,FALSE), each=2^pr[2])
    basis <- seq_len(2^n)[kp > 0] - 1

    rows <- vals + rep(basis,each=4)
    cols <- idx + rep(seq_len(2^(n-2))-1, each=4)
    vals <- c(1,-1,-1,1)

    M[cbind(rows, cols) + 1] <- vals

    ## get column names
    nms <- powerSet(seq_len(n)[-(pr+1)])
    nms <- sapply(nms, paste0, collapse=",")
    colnames(M)[idx+seq_len(2^(n-2))] <- paste(paste0(pr+1,collapse=","), "|", nms, sep="")

    ## advance column counter
    idx <- idx + length(basis)
  }

  return(M)
}

##' Test if an independence is represented in an imset
##'
##' @param imset an imset
##' @param ci a conditional independence
##' @param timeout a timeout for the linear program solver in seconds
##' @param sparse logical: use sparse matrix?
##'
##' @details Runs linear program suggested in Lindner (2012), Section 6.3.2. The
##' timeout variable defaults to 60 seconds: using 0 means there is no limit. If
##' the version of \code{lpSolve} 5.6.13.4.9000 is not installed, then the
##' timeout variable will not have any effect.
##'
##' @importFrom lpSolve lp
##'
##' @references Lindner (2012). \emph{Discrete optimisation in machine learning-learning of Bayesian network structures and conditional independence implication}, PhD Thesis, TUM.
##'
test_indep <- function (imset, ci, consMat, timeout=60L, sparse=FALSE) {

  n <- log2(length(imset))
  nc <- choose(n,2)*2^(n-2)
  nr <- 2^n

  v <- elem_imset(ci[[1]], ci[[2]], ci[[3]], n=n)

  if (missing(consMat)) {
    consMat <- gen_constraints(n, sparse=sparse)
    if (sparse) {
      consMat2 <- cbind(consMat@i+1, consMat@j+1, consMat@x)
      imset <- cbind(which(imset != 0), nc+1, -imset[imset != 0])
      diag <- cbind(nr + seq_len(nc+1), seq_len(nc+1), 1)
      consMat2 <- rbind(consMat2, imset, diag)
    }
    else {
      consMat2 <- cbind(matrix(consMat, ncol=nc), (-1)*imset)
      consMat2 <- rbind(consMat2, diag(nc+1))
    }
    consMat <- consMat2
  }
  # else sparse <- !is.matrix(consMat)
# for (i in seq_len(nr)) consMat2 <- cbind(consMat2, 0)

# col <- rep(seq_len(nc), times=diff(consMat@p))
# row <- consMat@i+1
# val <- c(1,-1,-1,1)
# spMat <- cbind(row,col,val)


## check if elem_imset is contained in structural imset
# which(u != 0)
#
# consMat2 <- rbind(consMat2, c(rep(1,nc),(-1)*u))
# consMat2 <- rbind(consMat2, do.call(cbind, c(rep(list(0), nc), list(cbind(diag(nr-1), -1)))))
# matrix(consMat2, nrow=nrow(consMat2))
#
# spMat <- rbind(spMat,
#                cbind(rep(nr+1,nc), seq_len(nc), 1),
#                cbind(rep(nr+1,nr), nc+seq_len(nr), (-1)*u))
#
# spMat <- rbind(spMat, cbind(nc+1+rep(seq_len(nr-1), 2),
#                             nc+rep(seq_len(nr-1),2),
#                             rep(c(1,-1),each=nr-1)))


  dirs <- c(rep("==", nr), rep(">=", nc+1))

  if (packageVersion("lpSolve") != '5.6.13.4.9000') {
    if (!missing(timeout)) message("Wrong version of lpSolve installed, so timeout will not work")

    if (sparse) {
      out <- lp(direction="min", rep(0,nc+1),
                dense.const = consMat,
                const.dir = dirs,
                const.rhs = c((-1)*v, rep(0,nc+1)),
                all.int = TRUE)
    }
    else {
      out <- lp(direction="min", rep(0,nc+1),
                const.mat = consMat,
                const.dir = dirs,
                const.rhs = c((-1)*v, rep(0,nc+1)),
                all.int = TRUE)
    }
  }
  else {
    if (sparse) {
      out <- lp(direction="min", rep(0,nc+1),
                dense.const = consMat,
                const.dir = dirs,
                const.rhs = c((-1)*v, rep(0,nc+1)),
                all.int = TRUE,
                timeout = timeout)
    }
    else {
      out <- lp(direction="min", rep(0,nc+1),
                const.mat = consMat,
                const.dir = dirs,
                const.rhs = c((-1)*v, rep(0,nc+1)),
                all.int = TRUE,
                timeout = timeout)
    }
  }
  if (out$status == 0) return(TRUE)
  else if (out$status == 2) return(FALSE)
  else if (out$status == 7) return(NA)
  else {
    message(paste0("Unknown status: ", out$status))
    return(NA)
  }
}

##' Check if the an imset is structural or combinatorial
##'
##' Uses linear programming to check whether imsets can be decomposed as
##' elementary imsets.
##'
##' @param imset an imset to be tested
##' @param timeout an integer giving the maximum run time for each linear program in seconds
##' @param sparse should sparse matrices be used?
##'
##' @details For the \code{timeout} operator to work we require version 5.6.13.4.9000
##' of \code{lpSolve} to be installed.
##'
##'
##' @export
is_combinatorial <- function (imset, timeout=60L, sparse=FALSE) {

  n <- log2(length(imset))
  nc <- choose(n,2)*2^(n-2)
  nr <- 2^n

  consMat <- gen_constraints(n, sparse=sparse)
  if (sparse) {
    consMat2 <- cbind(consMat@i+1, consMat@j+1, consMat@x)
    diag <- cbind(nr + seq_len(nc), seq_len(nc), 1)
    consMat2 <- rbind(consMat2, diag)
  }
  else {
    # consMat2 <- cbind(matrix(consMat, ncol=nc), (-1)*imset)
    consMat2 <- matrix(consMat, ncol=nc)
    consMat2 <- rbind(consMat2, diag(nc))
  }
  dirs <- c(rep("==", nr), rep(">=", nc))

  if (packageVersion("lpSolve") != '5.6.13.4.9000') {
    if (!missing(timeout)) message("Wrong version of lpSolve installed, so timeout will not work")

    if (sparse) {
      out <- lp(direction="min", rep(0,nc),
                dense.const = consMat2,
                const.dir = dirs,
                const.rhs = c(imset, rep(0,nc)),
                all.int = TRUE)
    }
    else {
      out <- lp(direction="min", rep(0,nc),
              const.mat = consMat2,
              const.dir = dirs,
              const.rhs = c(imset, rep(0,nc)),
              all.int = TRUE)
    }
  }
  else {
    if (sparse) {
      out <- lp(direction="min", rep(0,nc),
                dense.const = consMat2,
                const.dir = dirs,
                const.rhs = c(imset, rep(0,nc)),
                all.int = TRUE,
                timeout = timeout)
    }
    else {
      out <- lp(direction="min", rep(0,nc),
                const.mat = consMat2,
                const.dir = dirs,
                const.rhs = c(imset, rep(0,nc)),
                all.int = TRUE,
                timeout = timeout)
    }
  }

  if (out$status == 0) {
    ret <- TRUE
    sol <- out$solution
    names(sol) <- indep_labels(n)
    attr(ret, "indeps") <- sol[sol != 0]
    return(ret)
  }
  else if (out$status == 2) return(FALSE)
  else if (out$status == 7) return(NA)
  else {
    message(paste0("Unknown status: ", out$status))
    return(NA)
  }
}

indep_labels <- function(n) {
  if (n <= 1) return(character(0))
  len <- choose(n,2)*2^(n-2)
  if (len > 1e6) warning("Very long vector being created")

  prs <- combn(n, 2)
  out <- character(len)

  for (i in seq_len(choose(n,2))) {
    pr <- prs[,i]
    nms <- powerSet(seq_len(n)[-pr])
    nms <- sapply(nms, paste0, collapse=",")
    out[2^(n-2)*(i-1)+seq_len(2^(n-2))] <- paste(paste0(pr,collapse=","), ifelse(nchar(nms)>0,"|",""), nms, sep="")
  }

  out
}

##' @describeIn is_combinatorial test if imset is structural
##' @export
is_structural <- function (imset, timeout=60L, sparse=FALSE) {

  n <- log2(length(imset))
  nc <- choose(n,2)*2^(n-2)
  nr <- 2^n

  consMat <- gen_constraints(n, sparse=sparse)
  if (sparse) {
    consMat2 <- cbind(consMat@i+1, consMat@j+1, consMat@x)
    imset <- cbind(which(imset != 0), nc+1, -imset[imset != 0])
    diag <- cbind(nr + seq_len(nc+1), seq_len(nc+1), 1)
    consMat2 <- rbind(consMat2, imset, diag)
  }
  else {
    consMat2 <- cbind(matrix(consMat, ncol=nc), (-1)*imset)
    # consMat2 <- matrix(consMat, ncol=nc)
    consMat2 <- rbind(consMat2, diag(nc+1))
  }

  dirs <- c(rep("==", nr), rep(">=", nc), ">=")

  if (packageVersion("lpSolve") != '5.6.13.4.9000') {
    if (!missing(timeout)) message("Wrong version of lpSolve installed, so timeout will not work")

    if (sparse) {
      out <- lp(direction="min", c(rep(0,nc+1)),
                dense.const = consMat2,
                const.dir = dirs,
                const.rhs = c(rep(0,nr+nc),1),
                all.int = TRUE)
    }
    else {
      out <- lp(direction="min", c(rep(0,nc+1)),
                const.mat = consMat2,
                const.dir = dirs,
                const.rhs = c(rep(0,nr+nc),1),
                all.int = TRUE)
    }
  }
  else {
    if (sparse) {
      out <- lp(direction="min", c(rep(0,nc+1)),
                dense.const = consMat2,
                const.dir = dirs,
                const.rhs = c(rep(0,nr+nc),1),
                all.int = TRUE,
                timeout = timeout)
    }
    else {
      out <- lp(direction="min", c(rep(0,nc+1)),
                const.mat = consMat2,
                const.dir = dirs,
                const.rhs = c(rep(0,nr+nc),1),
                all.int = TRUE,
                timeout = timeout)
    }
  }

  if (out$status == 0) {
    ret <- TRUE

    sol <- out$solution[seq_len(nc)]
    names(sol) <- indep_labels(n)

    attr(ret, "indeps") <- sol[sol != 0]
    attr(ret, "k") <- last(out$solution)
    return(ret)
  }
  else if (out$status == 2) return(FALSE)
  else if (out$status == 7) return(NA)
  else {
    message(paste0("Unknown status: ", out$status))
    return(NA)
  }
}


##' Check if the standard imset for a graph defines the model
##'
##' @param graph an ADMG in the form of a \code{mixedgraph} object
##' @param u optionally, an imset to test the constraints from the graph on
##' (otherwise the 'standard' imset is used)
##' @param timeout an integer giving the maximum run time for each linear program in seconds
##' @param trace logical: should details be given of tests?
##' @param sparse logical: use sparse matrix?
##' @param rmv_edges logical: only use missing edges in the graph in the constraint matrix?
##'
##' @details The \code{timeout} argument requires version 5.6.13.4.9000 of \code{lpSolve} to be
##' installed.  If the function returns \code{NA} then no independence was found
##' to be not in the model but at least one of the linear programs timed out.
##'
##' @export
defines_mod <- function (graph, u, timeout=60L, trace=FALSE, sparse=FALSE,
                         rmv_edges=FALSE) {
  if (missing(u)) u <- standard_imset(graph)
  if (rmv_edges) {
    nind <- withEdgeList(morphEdges(graph, to="undirected"))$edges$undirected
  }
  else nind <- list()
  mod <- ADMGs2::localMarkovProperty(graph, split=TRUE)
  out <- TRUE

  if (packageVersion("lpSolve") != '5.6.13.4.9000' && !missing(timeout)) {
    message("Wrong version of lpSolve installed, so timeout will not work")
    timeout = NA
  }

  ## get constraint matrix
  n <- log2(length(u))
  consMat <- gen_constraints(n, sparse = sparse, reduce=nind)

  if (sparse) {
    nc <- ncol(consMat)
    nr <- nrow(consMat)

    consMat2 <- cbind(consMat@i+1, consMat@j+1, consMat@x)
    imset <- cbind(which(u != 0), nc+1, -u[u != 0])
    diag <- cbind(nr + seq_len(nc+1), seq_len(nc+1), 1)
    consMat <- rbind(consMat2, imset, diag)
  }
  else {
    consMat <- cbind(consMat, (-1)*u)
    consMat <- rbind(consMat, diag(ncol(consMat)))
  }

  ## now check each independence is represented
  for (k in seq_along(mod)) {
    if (trace) {
      cat("Testing independence: ")
      print(mod[[k]])
    }
    if (is.na(timeout)) ind_k <- test_indep(u, mod[[k]])
    else ind_k <- test_indep(u, mod[[k]], consMat=consMat, timeout=timeout, sparse=sparse)
    if (is.na(ind_k)) out <- NA
    else if(!ind_k) out <- FALSE

    if (isFALSE(out)) break
  }
  out
}

##' Check if the standard imset for a graph defines the model (deprecated)
##'
##' @param graph an ADMG in the form of a \code{mixedgraph} object
##' @param u optionally, an imset to test the constraints from the graph on
##' (otherwise the 'standard' imset is used)
##' @param timeout an integer giving the maximum run time for each linear program in seconds
##' @param trace logical: should details be given of tests?
##'
##' @details This requires version 5.6.13.4 or higher of \code{lpSolve} to be
##' installed.  If the function returns \code{NA} then no independence was found
##' to be not in the model but at least one of the linear programs timed out.
##'
##' This function has been superseded by \code{defines_mod}.
##'
##'
##' @export
definesMod <- function (graph, u, timeout=60L, trace=FALSE) {
  .Deprecated("defines_mod")

  if (missing(u)) u <- standard_imset(graph)
  mod <- ADMGs2::localMarkovProperty(graph, split=TRUE)
  out <- TRUE

  if (packageVersion("lpSolve") != '5.6.13.4.9000' && !missing(timeout)) {
    message("Wrong version of lpSolve installed, so timeout will not work")
    timeout = NA
  }

  for (k in seq_along(mod)) {
    if (trace) {
      cat("Testing independence: ")
      print(mod[[k]])
    }
    if (is.na(timeout)) ind_k <- test_indep(u, mod[[k]])
    else ind_k <- test_indep(u, mod[[k]], timeout=timeout)
    if (is.na(ind_k)) out <- NA
    else if(!ind_k) out <- FALSE

    if (isFALSE(out)) break
  }
  out
}

definesMod2 <- function (graph) {
  U <- standard_imset(graph)

  for (i in graph$v[-length(graph$v)]) for (j in graph$v[-seq_len(i)]) {
    if (j %in% adj(graph, i)) next

    non_m_seps <- list()
    for (k in powerSet(setdiff(graph$v, c(i,j)))) {
      # print(k)
      if (!m_sep(graph, i, j, k)) {
        non_m <- list(i,j,k)
        class(non_m) <- "ci"
        # print(non_m)
        non_m_seps <- c(non_m_seps, list(non_m))
      }
    }

    ## run linear programs
    for (k in seq_along(non_m_seps)) {
      print(non_m_seps[[k]])
      ind <- test_indep(U, non_m_seps[[k]])
      if (ind) stop()
    }
  }
  return(FALSE)
}

##' @describeIn is_combinatorial for a structural imset, get its degree
##' @export
imset_degree <- function(imset, timeout=60L) {
  out <- is_structural(imset, timeout=timeout)
  if (out) sum(attr(out, "indep"))
  else NA
}
