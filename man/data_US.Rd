\name{data_US}
\alias{data_US}
\docType{data}
\title{
Historical data for the US (GDP and Financial Conditions) from 1973:Q1 to 2022:Q3
}
\description{
data_euro contains:
- Quarterly annualized GDP, from 1973:Q1 to 2022:Q3
- National Financial Condition Index of the US, from 1973:Q1 to 2022:Q3
Sources :    https://www.chicagofed.org/research/data/nfci/current-data    https://fred.stlouisfed.org/series/A191RL1Q225SBEA
}
\usage{data("data_US")}
\format{
  A data frame with 200 observations on the following 4 variables.
  \describe{
    \item{\code{DATE}}{Vector of dates.}
    \item{\code{GDP}}{Vector of annualized PIB.}
    \item{\code{NFCI}}{Historical values of the National Financial Condition Index (NFCI).}
}
}
