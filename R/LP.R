gen_constraints <- function (n) {
  if (n < 0) stop("n must be non-negative")
  if (n <= 1) return(matrix(NA, nrow=0, ncol=n+1))

  prs <- combn(n, 2)
  M <- matrix(0, 2^n, choose(n,2)*2^(n-2))
  idx <- 0

  for (i in seq_len(choose(n,2))) {
    pr <- prs[,i]
    # rje::powerSet(seq_len(n-2))
    vals <- c(0,2^(pr[1]-1),2^(pr[2]-1),2^(pr[1]-1)+2^(pr[2]-1))

    basis <- seq_len(2^n)
    kp <- rep(c(TRUE,FALSE), each=2^(pr[1]-1))
    kp <- kp*rep(c(TRUE,FALSE), each=2^(pr[2]-1))
    basis <- basis[kp > 0]

    to_fill <- cbind(vals+rep(basis,each=4), rep(seq_along(basis), each=4))

    M[,idx+seq_along(basis)][to_fill] <- c(1,-1,-1,1)
    idx <- idx + length(basis)
  }

  M
}

##' Test if an independence is represented in an imset
##'
##' @param imset an imset
##' @param ci a conditional independence
##' @param timeout a timeout for the linear program solver in seconds
##'
##' @details Runs linear program suggested in Lindner (2012), Section 6.3.2. The
##' timeout variable defaults to 60 seconds: using 0 means there is no limit.
##'
##' @importFrom lpSolve lp
##'
##' @references Lindner (2012). \emph{Discrete optimisation in machine learning-learning of Bayesian network structures and conditional independence implication}, PhD Thesis, TUM.
##'
test_indep <- function (imset, ci, timeout=60L) {

  n <- log2(length(imset))
  nc <- choose(n,2)*2^(n-2)
  nr <- 2^n

  v <- elemImset(ci[[1]], ci[[2]], ci[[3]], n=n)

  conMat <- gen_constraints(n)
  conMat2 <- cbind(matrix(conMat, ncol=nc), (-1)*imset)
  conMat2 <- rbind(conMat2, diag(nc+1))


# for (i in seq_len(nr)) conMat2 <- cbind(conMat2, 0)

# col <- rep(seq_len(nc), times=diff(conMat@p))
# row <- conMat@i+1
# val <- c(1,-1,-1,1)
# spMat <- cbind(row,col,val)


## check if elemImset is contained in structural imset
# which(u != 0)
#
# conMat2 <- rbind(conMat2, c(rep(1,nc),(-1)*u))
# conMat2 <- rbind(conMat2, do.call(cbind, c(rep(list(0), nc), list(cbind(diag(nr-1), -1)))))
# matrix(conMat2, nrow=nrow(conMat2))
#
# spMat <- rbind(spMat,
#                cbind(rep(nr+1,nc), seq_len(nc), 1),
#                cbind(rep(nr+1,nr), nc+seq_len(nr), (-1)*u))
#
# spMat <- rbind(spMat, cbind(nc+1+rep(seq_len(nr-1), 2),
#                             nc+rep(seq_len(nr-1),2),
#                             rep(c(1,-1),each=nr-1)))


  dirs <- c(rep("==", nr), rep(">=", nc+1))

  out <- lp(direction="min", rep(0,nc+1),
            # dense.const = spMat,
            const.mat = conMat2,
            const.dir = dirs,
            const.rhs = c((-1)*v, rep(0,nc+1)),
            all.int = TRUE,
            timeout = timeout)
  if (out$status == 0) return(TRUE)
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
##' @param timeout an integer giving the maximum run time for each linear program in seconds
##'
##' @details This requires version 5.6.13.4 or higher of \code{lpSolve} to be
##' installed.  If the function returns \code{NA} then no independence was found
##' to be not in the model but at least one of the linear programs timed out.
##'
##' @export
definesMod <- function (graph, timeout=60L) {
  U <- standard_imset(graph)
  mod <- ADMGs2::localMarkovProperty(graph, split=TRUE)
  out <- TRUE
  for (k in seq_along(mod)) {

    ind_k <- test_indep(U, mod[[k]], timeout=timeout)
    if (is.na(ind_k)) out <- NA
    else if(!ind_k) out <- FALSE

    if (isFALSE(out)) break
  }
  out
}
