############################################################
# This is a dataset from Chan (2020, JBES) paper. 
# The detailed descriptions are given in Appendix B of Supplementary material
# available at: https://doi.org/10.1080/07350015.2018.1451336
############################################################
# List of variables:
# 
# Real gross domestic product
# Consumer price index
# Effective Federal funds rate
# M2 money stock
# Personal income
# Real personal consumption expenditure
# Industrial production index
# Civilian unemployment rate
# Housing starts
# Producer price index
# Personal consumption expenditures: chain-type price index
# Average hourly earnings: manufacturing
# MI money stock
# 10-Year Treasury constant maturity rate
# Real gross private domestic investment
# All employees: total nonfarm
# ISM manufacturing: PMI composite index
# ISM manufacturing: new orders index
# Business sector: real output per hour of all Persons
# Real stock prices (S& P 500 index divided by PCE index)

data_tmp  = read.csv("inst/varia/data/data_Q.csv")

us_macro_chan = ts(
  data_tmp[,c(1:3, 6:15, 17, 19:24)],
  start = c(1959, 4), 
  frequency = 4
)
colnames(us_macro_chan)  = c(
  "rgdp", 
  "cpi", 
  "FFR",
  "m2",
  "pinc",
  "rpce",
  "ip",
  "UR",
  "hs",
  "pci",
  "pce",
  "ahem",
  "mi",
  "TMR10Y",
  "rgpdi",
  "aetnf",
  "pmici",
  "noi",
  "bsro",
  "sp500"
)
save(us_macro_chan, file = "data/us_macro_chan.rda")
