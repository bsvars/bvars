# A 20-variable US macroeconomic system for the period 1959 Q4 – 2013 Q4

A system of 20 US macroeconomic aggregates used by Chan (2020).

## Usage

``` r
data(us_macro_chan)
```

## Format

A matrix and a `ts` object with time series of 217 observations on 20
variables:

- rgdp:

  Real gross domestic product

- cpi:

  Consumer price index

- FFR:

  Effective Federal funds rate

- m2:

  M2 money stock

- pinc:

  Personal income

- rpce:

  Real personal consumption expenditure

- ip:

  Industrial production index

- UR:

  Civilian unemployment rate

- hs:

  Housing starts

- pci:

  Producer price index

- pce:

  Personal consumption expenditures: chain-type price index

- ahem:

  Average hourly earnings: manufacturing

- mi:

  MI money stock

- TMR10Y:

  10-Year Treasury constant maturity rate

- rgpdi:

  Real gross private domestic investment

- aetnf:

  All employees: total nonfarm

- pmici:

  ISM manufacturing: PMI composite index

- noi:

  ISM manufacturing: new orders index

- bsro:

  Business sector: real output per hour of all Persons

- sp500:

  Real stock prices (S& P 500 index divided by PCE index

The series are used and described by Chan (2020) in Appendix B of
Supplementary Materials available at
\<doi:10.1080/07350015.2018.1451336\>.

## Source

FRED Economic Database, Federal Reserve Bank of St. Louis,
<https://fred.stlouisfed.org/>

## References

Chan (2020) Large Bayesian VARs: A Flexible Kronecker Error Covariance
Structure, Journal of Business and Economic Statistics, 38(1), 68–79,
\<doi:10.1080/07350015.2018.1451336\>.

## Examples

``` r
data(us_macro_chan)   # upload the data
```
