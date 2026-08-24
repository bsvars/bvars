
set.seed(1)
run_no1     <- rmatnorm1(matrix(0, 2, 3), diag(2), diag(3))

set.seed(1)
run_no2     <- rmatnorm1(matrix(0, 2, 3), diag(2), diag(3))

expect_true(
  all(dim(run_no1) == dim(run_no2)),
  info = "rmatnorm1: same output dimentions for normal and pipe workflow."
)

expect_identical(
  run_no1[1,1], run_no2[1,1],
  info = "rmatnorm1: identical for normal and pipe workflow."
)

expect_error(
  rmatnorm1(matrix(0, 2, 2), diag(2), diag(3)),
  info = "rmatnorm1: dimensions don't match."
)

expect_error(
  rmatnorm1(matrix(0, 2, 3), diag(3), diag(3)),
  info = "rmatnorm1: dimensions don't match."
)

expect_error(
  rmatnorm1(matrix(0, 2, 3), diag(2), diag(2)),
  info = "rmatnorm1: dimensions don't match."
)

M <- matrix(0, 2, 3)
M[1,1] <- NA
expect_error(
  rmatnorm1(M, diag(2), diag(3)),
  info = "rmatnorm1: misspecified argument."
)

