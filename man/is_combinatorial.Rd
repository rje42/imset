% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combinatorial.R
\name{is_combinatorial}
\alias{is_combinatorial}
\title{Recursive function to check if an imset if combinatorial}
\usage{
is_combinatorial(imset, use_char = TRUE, trace = FALSE)
}
\arguments{
\item{imset}{imset to be checked}

\item{use_char}{logical: should characteristic imset method be used?}

\item{trace}{logical: show detail}
}
\description{
Recursive function to check if an imset if combinatorial
}
\details{
This function recursively tries to remove conditional independences
by greedy searching.  It looks at the largest set still present, and then
tries all pairs still in the characteristic imset that are contained in
this set as variables to be made conditionally independent.  The associated
elementary characteristic imset is then subtracted, until we either reach
the zero imset (in which case the imset is combinatorial), or we see a
negative value.  In the latter case, we reject that sequence and try
another branch.  Once all branches are exhausted the algorithm fails.
}
\seealso{
\code{\link{defines_mod}}
}