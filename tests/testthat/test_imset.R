elem <- elem_imset(1,2,3,n=4)
c_elem <- char_imset(elem)
elem2 <- standard_imset(c_elem)

test_that("imsets work", {
  expect_equal(elem, elem2)
})

gr5 <- makeGraphCycle(5, "bidirected")
ch5 <- char_imset(gr5)
st5 <- standard_imset(gr5)

standard_imset(ch5)
test_that("graphical imsets work", {
  expect_equal(char_imset(st5), ch5)
  expect_equal(standard_imset(ch5), st5)
})
