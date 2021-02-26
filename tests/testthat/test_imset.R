elem <- elemImset(1,2,3,n=4)
c_elem <- char_imset(elem)
elem2 <- standard_imset(c_elem)

test_that("imsets work", {
  expect_equal(elem, elem2)
})
