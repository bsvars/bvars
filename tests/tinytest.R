
if (requireNamespace("tinytest", quietly = TRUE)) {
  home <- length(unclass(packageVersion("bvars"))[[1]]) == 4
  tinytest::test_package("bvars", at_home = home)
}

