##' Counterexample from Hemmecke et al. (2008)
##'
##' This is an imset that is structural but not combinatorial.
##' Code to recreate it is given below.
##'
##' @examples
##' is_combinatorial(hem08)
##' is_structural(hem08)
##'
##' ## code to recreate hem08
##' u <- elem_imset(4, 5, 2) + elem_imset(4, 5, 3) +
##'   elem_imset(1, 3, 4) + elem_imset(1, 2, 5) +
##'   elem_imset(2, 5, c(1,4)) + elem_imset(4, 3, c(1,5)) +
##'   elem_imset(1, 4, c(2,3)) + elem_imset(1, 5, c(2,3)) +
##'   elem_imset(1, 4, c(3)) + elem_imset(1, 5, c(2)) +
##'   elem_imset(2, 3, 4) + elem_imset(2, 3, 5) +
##'   elem_imset(4, 3, 1:2) + elem_imset(2, 5, c(1,3)) +
##'   elem_imset(1, 3, 4:5) + elem_imset(2, 1, 4:5)
##'
##' @references Hemmecke et al. (2008) Three Counter-Examples on Semi-Graphoids, \emph{Combinatorics Probability and Computing}, 17(2):239-257.
"hem08"

##' Pair of bidirected graphs
##'
##' Two graphs, \code{gr1a} and \code{gr1b}, that differ by one edge, but one
##' (\code{gr1a}) is perfectly Markovian with respect to its standard imset,
##' and the other (\code{gr1b}) is not.
##'
##' @aliases gr1a gr1b
##'
"gr1a"

##' Bidirected graphs on 5 and 6 vertices
##'
##' List of all bidirected connected graphs on 5 and 6 unlabelled vertices.
##'
##' @aliases bidi5 bidi6
"bidi5"

##' MAGs that do not defined the model over 5 and 6 vertices
##'
"fail56"
