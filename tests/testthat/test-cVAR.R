library(ggplot2)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)

#---- SW2001 ----

## prepare data
rm(list=ls())
load_all()
data(SW2001)
data.1 <- SW2001
head(data.1)

## estimate VAR, identify and compute IRF
myvar <- cVAR$new(data = data.1, p.lag = 4)
head(myvar$data)
myvar$est()
myvar$identify(method = "recursive", boot.method = "block")
IRF <- myvar$tool("IRF", 24, 0.90)

## plot preparation
shock.names <- c("unemployment", "inflation", "interest rate")
var.names <- c("unemp", "infl", "FFR")
xaxis <- 0:24

## base R plot
par(mfrow = c(3, 3))
par(family="serif", cex.main = 1.5, cex.lab = 1.2, cex.axis=1.2)
for (i in 1:3) {
  for (j in 1:3) {
    shk.name <- shock.names[j]
    var.name <- var.names[i]
    base <- IRF$base[i, j, ]
    ub <- IRF$ub[i, j, ]
    lb <- IRF$lb[i, j, ]
    y.range <- range(c(base, lb, ub))
    plot(xaxis, base,
         type = "l", col = "blue", ylim = y.range, lwd=3,
         xlab = "Months", ylab = "Percent", main = paste(var.name,"to",shk.name,"shock")
    )

    polygon(c(xaxis, rev(xaxis)), c(ub, rev(lb)), col = "#C0C0C0AE", border = NA)
    lines(xaxis, base, col = "blue", lty = 1, lwd = 3)
    abline(h = 0, lwd = 3, lty = 2)
  }
}

## ggplot2
plot_data <- expand.grid(
  Shock = shock.names,
  Variable = var.names,
  Month = xaxis
) %>%
  mutate(
    Base = as.vector(IRF$base),
    Upper = as.vector(IRF$ub),
    Lower = as.vector(IRF$lb)
  )

# Create the plot
ggplot(plot_data, aes(x = Month, y = Base)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#C0C0C0AE", alpha = 0.5) +  # Shaded confidence band
  geom_line(color = "blue", size = 1) +  # Line for the base
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +  # Horizontal line at 0
  facet_grid(Variable ~ Shock, scales = "free", space = "free") +  # Multi-panel layout
  labs(
    x = "Months",
    y = "Percent",
    title = "Impulse Response Functions",
    caption = "solid line is base IRF, shaded area is 90% block bootstrap CI"
  ) +
  theme_classic(base_family = "serif", base_size = 14) +
  theme(
    strip.text = element_text(size = 12),  # Panel titles
    plot.title = element_text(size = 16, hjust = 0.5)
  )

# ---- GK 2015 ----

load_all()
data("GK2015")
data.2 <- GK2015
head(GK2015)
iv <- as.matrix(data.2$ff4_tc)
colnames(iv) <- c("gs1")
variv <- cVAR$new(data = data.2[,1:5], p.lag = 12, data.IV = iv)
variv$est()
variv$identify(method="IV", boot.method = "block")
IRF <- variv$tool("IRF", 48, 0.90)

## plot preparation
shk.name <- "one year rate"
var.names <- c("CPI", "Industry Production", "one year rate", "Excess bond premium")
xaxis <- 0:48

## base R plot
par(mfrow = c(4, 1))
par(family="serif", cex.main = 1.5, cex.lab = 1.2, cex.axis=1.2)

for (j in 1:4) {
  var.name <- var.names[j]
  base <- IRF$base[j, 3, ]*0.2
  ub <- IRF$ub[j, 3, ]*0.2
  lb <- IRF$lb[j, 3, ]*0.2
  y.range <- range(c(base, lb, ub))
  plot(xaxis, base,
       type = "l", col = "blue", ylim = y.range, lwd=3,
       xlab = "Months", ylab = "Percent", main = paste(var.name,"to",shk.name,"shock")
  )
  polygon(c(xaxis, rev(xaxis)), c(ub, rev(lb)), col = "#C0C0C0AE", border = NA)
  lines(xaxis, base, col = "blue", lty = 1, lwd = 3)
  abline(h = 0, lwd = 3, lty = 2)
}
