entry_names.numeric <- function(x, vnames, sep=",") {
  entry_names.list(posToSubset(x), vnames, sep)
}

entry_names.list <- function(lst, vnames, sep=",") {
  sapply(lst, function(x) ifelse(length(x) > 0, paste(vnames[x], collapse=sep), "()"))
}
