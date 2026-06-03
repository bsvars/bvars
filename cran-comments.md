## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.
* performed all checks using `usethis::use_release_issue()` and we pass
* performed all checks using `devtools::check(remote = TRUE, manual = TRUE)` and we pass except for:
`Package in Depends/Imports which should probably only be in LinkingTo: ‘RcppArmadillo’`
We investigated and the current setting is the only one that passes all other checks.
* all checks from `devtools::check_win_devel()` done!
* GitHub repo actions checking the package all pass.