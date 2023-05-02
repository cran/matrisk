\name{data_euro}
\alias{data_euro}
\docType{data}
\title{
Historical data for the eurozone (GDP and Financial Conditions) from 2008:Q4 to 2022:Q3
}
\description{
data_euro contains:
- Quarterly annualized GDP, from 2008:Q4 to 2022:Q3
- Financial Condition Index of the euro Area, from 2008:Q4 to 2022:Q3
- Composite Indicator of Systemic Stress, from 2008:Q4 to 2022:Q3
Sources :    https://sdw.ecb.europa.eu/browseExplanation.do?node=9689686    https://webstat.banque-france.fr/ws_wsen/browseSelection.do?node=DATASETS_FCI    https://fred.stlouisfed.org/series/CLVMEURSCAB1GQEA19
}
\usage{data("data_euro")}
\format{
  A data frame with 57 observations on the following 4 variables.
  \describe{
    \item{\code{DATE}}{Vector of dates.}
    \item{\code{GDP}}{Vector of annualized PIB.}
    \item{\code{FCI}}{Historical values of the Financial Condition Index (FCI).}
    \item{\code{CISS}}{Historical values of the Composite Indicator of Systemic Stress (CISS).}
}
}
