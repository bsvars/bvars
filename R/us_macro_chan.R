
#' @title A 20-variable US macroeconomic system for the period 1959 Q4 -- 2013 Q4
#'
#' @description A system of 20 US macroeconomic aggregates used by Chan (2020).
#'
#' @usage data(us_macro_chan)
#' 
#' @format A matrix and a \code{ts} object with time series of 217 observations 
#' on 20 variables:
#' \describe{
#'  \item{rgdp}{Real gross domestic product}
#'  \item{cpi}{Consumer price index}
#'  \item{FFR}{Effective Federal funds rate}
#'  \item{m2}{M2 money stock}
#'  \item{pinc}{Personal income}
#'  \item{rpce}{Real personal consumption expenditure}
#'  \item{ip}{Industrial production index}
#'  \item{UR}{Civilian unemployment rate}
#'  \item{hs}{Housing starts}
#'  \item{pci}{Producer price index}
#'  \item{pce}{Personal consumption expenditures: chain-type price index}
#'  \item{ahem}{Average hourly earnings: manufacturing}
#'  \item{mi}{MI money stock}
#'  \item{TMR10Y}{10-Year Treasury constant maturity rate}
#'  \item{rgpdi}{Real gross private domestic investment}
#'  \item{aetnf}{All employees: total nonfarm}
#'  \item{pmici}{ISM manufacturing: PMI composite index}
#'  \item{noi}{ISM manufacturing: new orders index}
#'  \item{bsro}{Business sector: real output per hour of all Persons}
#'  \item{sp500}{Real stock prices (S& P 500 index divided by PCE index}
#' }
#' 
#' The series are used and described by Chan (2020) in Appendix B of 
#' Supplementary Materials available at <doi:10.1080/07350015.2018.1451336>.
#' 
#' @references 
#' Chan (2020) Large Bayesian VARs: A Flexible Kronecker Error Covariance Structure,
#' Journal of Business and Economic Statistics, 38(1), 68--79,
#' <doi:10.1080/07350015.2018.1451336>.
#' 
#' @source 
#' FRED Economic Database, Federal Reserve Bank of St. Louis, 
#' \url{https://fred.stlouisfed.org/}
#' 
#' @examples 
#' data(us_macro_chan)   # upload the data
"us_macro_chan"
