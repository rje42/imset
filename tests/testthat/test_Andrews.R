data(fail56)
out <- logical(56)

for (i in seq_len(56)) {
  u <- standard_imset(1 - imset:::NIE(fail56[[i]])$inc)
  out[i] <- defines_mod(fail56[[i]], u=u, sparse = TRUE, rmv_edges = TRUE)
}

testthat::test_that("Inclusion imset defines model", {
  expect_equal(out, rep(TRUE, 56))
})
