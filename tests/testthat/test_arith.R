elem_4 <- elem_imset(1,2,3,n=4)
elem_3 <- elem_imset(1,2,3)
elem_4a <- elem_imset(3,4)


test_that("imsets work", {
  expect_equal(elem_4a - elem_4, elem_4a - elem_3)
  expect_equal(elem_3 - elem_4, 0*elem_4)
})
