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

#' Uhlig (2005) dataset
#'
#' The data used in Uhlig (2005) paper *What are the effects of monetary policy on output? Results from an agnostic identification procedure*
#'
#' @format ## Uhlig2005
#'
#' The data used in Uhlig (2005) paper *What are the effects of monetary policy on output? Results from an agnostic identification procedure*, monthly data, 515 obs of 7 variables
#' \describe{
#'    \item{Date}{month indicator, 1965m1 - 2007m11}
#'    \item{GDPC1}{real GDP of USA}
#'    \item{GDPDEF}{GDP deflator}
#'    \item{CPRINDEX}{CPI}
#'    \item{TOTRESNS}{total reserves}
#'    \item{BOGNONBR}{non borrowed reserves}
#'    \item{FEDFUNDS}{federal funds rate}
#' }
#' @usage data(Uhlig2005)
#' @source https://doi.org/10.1016/j.jmoneco.2004.05.007
"Uhlig2005"

#' Kilian and Murphy (2012) dataset
#'
#' The data used in Kilian and Murphy (2012) paper *Why agnostic sign restrictions are not enough: understanding the dynamics of oil market var models*
#'
#' @format ## Kilian2012
#' The data used in Kilian and Murphy (2012) paper *Why agnostic sign restrictions are not enough: understanding the dynamics of oil market var models*, monthly data, 540 obs of 4 variables
#' \describe{
#'    \item{Date}{month indicator, 1971m1 - 2015m12}
#'    \item{Oil.Production.Growth}{global oil production growth}
#'    \item{Economic.Activity.Index}{global economic activity index}
#'    \item{Real.Oil.Price}{global real oil price}
#' }
#' @usage data(Kilian2012)
#' @source  https://doi.org/10.1111/j.1542-4774.2012.01080.x
"Kilian2012"

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
#'    \item{gs1}{1 year USA government bond rate, the monetary policy variable}
#'    \item{ebp}{excess bond premium}
#'    \item{ff4_tc}{three month ahead monthly fed funds futures, the IV}
#' }
#' @usage data(GK2015)
#' @source http://doi.org/10.3886/E114082V1
"GK2015"

