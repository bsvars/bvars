
#' @title Samples random numbers from the inverted gamma 2 distribution
#' @param n a positive integer number of draws
#' @param s a positive scalar scale parameter
#' @param nu a positive integer shape parameter
#' @return  a vector with \code{n} independent draws from the inverted gamma 2 distribution
#' @export
rig2 <- function(n, s, nu) {
    stopifnot(n > 0 & s > 0 & nu > 0)
    stopifnot(n %% 1 == 0)
    qq = .Call(`_bvars_rig2`, n, s, nu)
    return(qq)
}

