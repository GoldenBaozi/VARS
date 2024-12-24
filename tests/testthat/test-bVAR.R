library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(ggplot2)
library(microbenchmark)
# ---- Kilian (2012) ----
rm(list=ls())
load_all()
data("Kilian2012")
head(Kilian2012)
colnames(Kilian2012) <- c("Date", "Oil Production Growth", "Economic Activity Index", "Real Oil Price")
oilvar <- bVAR$new(data = Kilian2012, p.lag = 24)
oilvar$est(method = "conjugate")

SR <- list(
  list("Oil Production Growth", c("Oil Production Growth", "Economic Activity Index"), 0, 1),
  list("Oil Production Growth", "Real Oil Price", 0, -1),
  # list("Oil Production Growth", "Oil Production Growth", 1, 1),
  list("Economic Activity Index", c("Oil Production Growth", "Economic Activity Index", "Real Oil Price"), 0, 1),
  list("Real Oil Price", c("Oil Production Growth", "Real Oil Price"), 0, 1),
  list("Real Oil Price", "Economic Activity Index", 0, -1)
)

EBR <- list(
  list("Economic Activity Index", c("Oil Production Growth", "Real Oil Price"), 0, 0.0258, NA),
  list("Real Oil Price", c("Oil Production Growth", "Real Oil Price"), 0, 0.0258, NA)
)

NSR <- list(
  # list("Oil Production Growth", "contribution", as.Date("1990-08-01"), as.Date("1990-08-01"), "Oil Production Growth", -1, "most"),
  list("Economic Activity Index", "contribution", as.Date("1990-08-01"), as.Date("1990-08-01"), "Real Oil Price", -1, "least")
)
str(oilvar$.__enclos_env__$private$get.restrictions(SR, EBR, NSR))

oilvar$identify.seq(SR = SR, EBR = EBR, NSR = NSR, draw = 100, M = 1000)

IRF.1 <- oilvar$IRF.compute.batch(18, 0.68)
IRF.2 <- oilvar$IRF.compute.NSR(18, 0.68)
var.names <- c("Oil production", "Economic Activity Index", "Real Oil Price")
shock.names <- c("Oil Supply", "Aggregate Demand", "Oil-specific Demand")
xaxis <- 0:18

png("./IRF_oil.png", width=12, height=9, units="in", res=300)
par(mfrow = c(3, 3))
par(family="serif", cex.main = 1.5, cex.lab = 1.2, cex.axis=1.2)
for (i in 1:3) {
  for (j in 1:3) {
    shk.name <- shock.names[i]
    var.name <- var.names[j]
    med <- IRF.1$med[j, i, ]
    ub <- IRF.1$ub[j, i, ]
    lb <- IRF.1$lb[j, i, ]
    med.1 <- IRF.2$avg[j, i, ]
    ub.1 <- IRF.2$ub[j, i, ]
    lb.1 <- IRF.2$lb[j, i, ]
    if (i == 1) {
      med = -med
      ub = -ub
      lb = -lb
      med.1 = -med.1
      ub.1 = -ub.1
      lb.1 = -lb.1
    }
    y.range <- range(c(med, med.1, lb, lb.1, ub, ub.1))
    plot(xaxis, med,
         type = "l", col = "blue", ylim = y.range, lwd=3,
         xlab = "Months", ylab = "Percent", main = paste(var.name,"to",shk.name,"shock")
    )

    polygon(c(xaxis, rev(xaxis)), c(ub, rev(lb)), col = "#C0C0C0AE", border = NA)
    # lines(xaxis, median, col = "blue", lty = 1, lwd = 3)
    polygon(c(xaxis, rev(xaxis)), c(ub.1, rev(lb.1)), col = "#FFC0CBAB", border = NA)
    lines(xaxis, med.1, col = "red", lty = 1, lwd = 3)
    lines(xaxis, med, col = "blue", lty = 1, lwd = 3)
    abline(h = 0, lwd = 3, lty = 2)
  }
}
dev.off()

## ---- numerical test ----
oilvar$get.result.tmp <- FALSE

test_sample <- function(varobj) {
  varobj$identify.seq(SR = SR, EBR = EBR, draw = 100, M = 1000)
  varobj$get.result.tmp <- FALSE
}

mpvar$get.result.tmp <- TRUE

test_NSR <- function(varobj) {
  varobj$identify.seq(SR = SR, EBR = EBR, NSR = NSR, draw = 100, M = 1000)
}


microbenchmark(test_sample(oilvar), times = 5)

oilvar$get.result.tmp <- TRUE
microbenchmark(test_NSR(oilvar), times = 5)

# ---- Uhlig (2005) ----
data("Uhlig2005")
head(Uhlig2005)
mpvar <- bVAR$new(data = Uhlig2005, p.lag = 12, cons = 0)
mpvar$est(method = "conjugate")
SR <- list(
  list("FEDFUNDS", c("FEDFUNDS"), 0, 1),
  list("FEDFUNDS", c("GDPDEF", "CPRINDEX", "BOGNONBR"), 0, -1),
  list("FEDFUNDS", c("FEDFUNDS"), 1, 1),
  list("FEDFUNDS", c("GDPDEF", "CPRINDEX", "BOGNONBR"), 1, -1),
  list("FEDFUNDS", c("FEDFUNDS"), 2, 1),
  list("FEDFUNDS", c("GDPDEF", "CPRINDEX", "BOGNONBR"), 2, -1),
  list("FEDFUNDS", c("FEDFUNDS"), 3, 1),
  list("FEDFUNDS", c("GDPDEF", "CPRINDEX", "BOGNONBR"), 3, -1),
  list("FEDFUNDS", c("FEDFUNDS"), 4, 1),
  list("FEDFUNDS", c("GDPDEF", "CPRINDEX", "BOGNONBR"), 4, -1),
  list("FEDFUNDS", c("FEDFUNDS"), 5, 1),
  list("FEDFUNDS", c("GDPDEF", "CPRINDEX", "BOGNONBR"), 5, -1)
)

NSR <- list(
  list("FEDFUNDS", "sign", as.Date("1979-10-01"), as.Date("1979-10-01"), 1),
  list("FEDFUNDS", "contribution", as.Date("1979-10-01"), as.Date("1979-10-01"), "FEDFUNDS", 1, "overwhelm")
)

mpvar$identify.seq(SR = SR, NSR = NSR, draw = 1000, M = 1000)

IRF.1 <- mpvar$IRF.compute.batch(60, 0.68)
IRF.2 <- mpvar$IRF.compute.NSR(60, 0.68)
var.names <- c("GDP", "GDP deflator", "Commodity Price", "Total Reserve", "Non-Borrowed Reserve", "Fed Funds Rate")
xaxis <- 0:60

png("./monetary.png", width=12, height=6, units = "in", res=300)
par(mfrow = c(2, 3))
par(family="serif", cex.main = 1.5, cex.lab = 1.2, cex.axis=1.2)
for (i in 1:6) {
  j = 6
  shk.name <- "MP"
  var.name <- var.names[i]
  med <- IRF.1$med[i, j, ]*125
  ub <- IRF.1$ub[i, j, ]*125
  lb <- IRF.1$lb[i, j, ]*125
  med.1 <- IRF.2$avg[i, j, ]*71.6
  ub.1 <- IRF.2$ub[i, j, ]*71.6
  lb.1 <- IRF.2$lb[i, j, ]*71.6
  y.range <- range(c(med, med.1, lb, lb.1, ub, ub.1))
  plot(xaxis, med,
       type = "l", col = "blue", ylim = y.range, lwd=3,
       xlab = "Months", ylab = "Percent", main = paste(var.name,"to",shk.name,"shock")
  )

  polygon(c(xaxis, rev(xaxis)), c(ub, rev(lb)), col = "#C0C0C0AE", border = NA)
  # lines(xaxis, median, col = "blue", lty = 1, lwd = 3)
  polygon(c(xaxis, rev(xaxis)), c(ub.1, rev(lb.1)), col = "#FFC0CBAB", border = NA)
  lines(xaxis, med.1, col = "red", lty = 1, lwd = 3)
  lines(xaxis, med, col = "blue", lty = 1, lwd = 3)
  abline(h = 0, lwd = 3, lty = 2)
}
dev.off()

## ---- numerical test ----
mpvar$get.result.tmp <- FALSE

test_sample <- function(varobj) {
  varobj$identify.seq(SR = SR, draw = 100, M = 1000)
  varobj$get.result.tmp <- FALSE
}

microbenchmark(test_sample(mpvar), times = 10)

mpvar$get.result.tmp <- TRUE

test_NSR <- function(varobj) {
  varobj$identify.seq(SR = SR, NSR = NSR, draw = 100, M = 1000)
}

microbenchmark(test_NSR(mpvar), times = 10)
