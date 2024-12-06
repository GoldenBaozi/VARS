#' Stock and Watson (2001) dataset
#'
#' The data used in Stock and Watson (2001) paper *Vector Autoregressions*.
#'
#' @format ## SW2001
#' The data used in Stock and Watson (2001) paper *Vector Autoregressions*, quarterly data, 164 obs. of 4 variables.
#' \describe{
#'    \item{date}{quarter indicator, 1960q1 - 2000q4}
#'    \item{unemp}{USA unemployment rate}
#'    \item{infl}{USA inflation rate}
#'    \item{ffr}{USA federal funds rate}
#' }
#' @usage data(SW2001)
#' @source https://www.aeaweb.org/articles?id=10.1257/jep.15.4.101
"SW2001"

#' Gertler and Karadi (2015) dataset
#'
#' The data used in Gertler and Karadi (2015) paper *Monetary Policy Surprises, Credit Costs, and Economic Activity*.
#'
#' @format ## GK2015
#' The data used in Gertler and Karadi (2015) paper *Monetary Policy Surprises, Credit Costs, and Economic Activity*, monthly data, 396 obs of 6 variables.
#' \describe{
#'    \item{date}{month indicator, 1979m7 - 2012m6}
#'    \item{logcpi}{log of USA CPI}
#'    \item{logip}{log of USA industrial production}
#'    \item{gs1}{1 year USA government bond rate, the policy variable}
#'    \item{ebp}{excess bond premium}
#'    \item{ff4_tc}{three month ahead monthly fed funds futures, the IV}
#' }
#' @usage data(GK2015)
#' @source http://doi.org/10.3886/E114082V1
"GK2015"
