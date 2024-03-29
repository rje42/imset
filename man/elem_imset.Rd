% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elementary_imsets.R, R/utility.R
\name{elem_imset}
\alias{elem_imset}
\alias{elem_imset.default}
\alias{elem_imset.ci}
\alias{sei2ci}
\alias{ci2sei}
\title{Get a (semi-)elementary imset}
\usage{
elem_imset(...)

elem_imset.default(A, B, C = integer(0), n = max(c(A, B, C)), check = TRUE)

elem_imset.ci(ci, n)

sei2ci(imset)

ci2sei(ci, n = max(vars))
}
\arguments{
\item{A, B, C}{disjoint subsets of 1,..,n}

\item{n}{number of variables involved}

\item{check}{logical: check entries are valid?}

\item{ci}{object of class \code{ci}}

\item{imset}{semi-elementary imset}
}
\description{
Returns a semi-elementary imset \eqn{u_{\langle A, B | C \rangle}}{u<A, B | C>}
(as defined in Studeny, 2005) given numeric arguments or a \code{ci} object.
}
\details{
Returns an imset, a vector of length \eqn{2^n}
with all but four entries zero.
}
\section{Functions}{
\itemize{
\item \code{elem_imset.default()}: Method for integer vectors

\item \code{elem_imset.ci()}: Method for \code{ci} object

\item \code{sei2ci()}: 

\item \code{ci2sei()}: convert conditional independence to semi-elementary imset

}}
