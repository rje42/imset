u <- elem_imset(1,2) + elem_imset(1,4,2) + elem_imset(1,3)

test_that("is_combinatorial works", {
  expect_equivalent(is_combinatorial(u, sparse=FALSE), TRUE)
  expect_equivalent(is_combinatorial(u, sparse=TRUE), TRUE)
  expect_equivalent(is_combinatorial(u - elem_imset(2,4), sparse=FALSE), FALSE)
  expect_equivalent(is_combinatorial(u - elem_imset(2,4), sparse=TRUE), FALSE)
})

g5ch <- makeGraphChain(5, "b")
g5cy <- makeGraphCycle(5, "b")

test_that("defines_model works", {
  expect_equivalent(defines_mod(g5ch, sparse=FALSE), TRUE)
  expect_equivalent(defines_mod(g5ch, sparse=TRUE), TRUE)
  expect_equivalent(defines_mod(g5cy, sparse=FALSE), FALSE)
  expect_equivalent(defines_mod(g5cy, sparse=TRUE), FALSE)
  })
