
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bvar$new(us_macro_chan)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
fv                  <- compute_shocks(run_no1)

set.seed(1)
suppressMessages(
  fv2               <- us_macro_chan |>
    specify_bvar$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    compute_shocks()
)



expect_true(
  all(dim(fv) == dim(fv2)),
  info = "compute_shocks: same output dimentions for normal and pipe workflow."
)

expect_identical(
  fv[1,1,1], fv2[1,1,1],
  info = "compute_shocks: identical for normal and pipe workflow."
)
