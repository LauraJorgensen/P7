library(tidyverse)
library(Matrix)
library(ggplot2)
library(SQUAREM)

source("Desktop/Uni Mat/NMF - R code/minVolNMFfit.R")
source("Desktop/Uni Mat/NMF - R code/NMFfit.R")
source("Desktop/Uni Mat/NMF - R code/SFS-master/R/NMFPois.R")

invisible(sapply(list.files("/Users/amaliebechpoulsen/Desktop/Uni Mat/NMF - R code/SFS-master/R", full.names = TRUE, pattern = "\\.R$"), source))
Rcpp::sourceCpp("/Users/amaliebechpoulsen/Desktop/Uni Mat/NMF - R code/SFS-master/src/sampleSFS.cpp")

obj_name <- load("Desktop/Uni Mat/NMF - R code/LACA24.Rdata")

LACA <- as_tibble(get("V"))

M = as.matrix(LACA[, order(colSums(LACA))])

testk <- k.choice(trials = 5:15, V = M)

testk$bicall

df <- data.frame(k = as.integer(names(testk$bicall)), bic = as.numeric(testk$bicall))

ggplot(df, aes(x = k, y = bic)) +
  geom_line(linewidth = 0.8, colour = "dodgerblue4") +
  geom_point(size = 1.5) + geom_point(data = df[df$k == testk$k, ], size = 3, shape = 1, colour = "black") +
  scale_x_continuous(breaks = df$k) +
  labs(x = expression(italic(k)), y = "BIC")+
  theme_bw() +
  theme(panel.grid.minor = element_blank())
