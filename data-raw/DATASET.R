library(R.matlab)
SW2001 <- read.csv("./inst/extdata/SW2001_Data.csv")
head(SW2001)
usethis::use_data(SW2001)

GK2015 <- read.csv("./inst/extdata/GK2015_Data.csv")
head(GK2015)
use_data(GK2015)

Kilian2012 <- read.csv("./inst/extdata/Kilian_Data_Updated.csv")
head(Kilian2012)
str(Kilian2012)
Kilian2012 <- Kilian2012[, 2:5]
head(Kilian2012)
use_data(Kilian2012)

df <- readMat("./inst/extdata/Uhlig_Data_Updated.mat")
head(df)
names <- vector("character", 6)
for (i in 1:6) {
  names[i] <- df$varNames[1, i][[1]][[1]][1, 1]
}
dates <- as.Date(df$dates - 719529)
Uhlig2005 <- data.frame(
  Date = dates,
  df$data
)
colnames(Uhlig2005)[2:7] <- names
head(Uhlig2005)
use_data(Uhlig2005)
