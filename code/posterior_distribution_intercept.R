library(tidyverse)
library(BEST) # posterior distribution plot

post <- posterior_samples(bmod1)

plotPost(
    post$b_intercept, showCurve = TRUE, showMode = TRUE,
    credMass = NULL, xlab = "", xlim = c(0, 0.6)
    )

intercept <- post$b_intercept
mu <- mean(intercept)
sigma <- sd(intercept)

# quantiles_intercept <- quantile(
#     intercept, probs = c(0.05, 0.125, 0.25, 0.75, 0.875, 0.95)
#     )

hdi_50 <- hdi(intercept, credMass = 0.5) %>% as.numeric
hdi_80 <- hdi(intercept, credMass = 0.8) %>% as.numeric
hdi_95 <- hdi(intercept, credMass = 0.95) %>% as.numeric

# abline(v = quantiles_intercept[1], lty = 2)
# abline(v = quantiles_intercept[2], lty = 2)
# abline(v = quantiles_intercept[3], lty = 2)
# abline(v = quantiles_intercept[4], lty = 2)
# abline(v = quantiles_intercept[5], lty = 2)
# abline(v = quantiles_intercept[6], lty = 2)

abline(v = hdi_50[1], lty = 2)
abline(v = hdi_50[2], lty = 2)
abline(v = hdi_80[1], lty = 2)
abline(v = hdi_80[2], lty = 2)
abline(v = hdi_95[1], lty = 2)
abline(v = hdi_95[2], lty = 2)

library(export) # export to ppt
graph2ppt(file = "plotpost.pptx", width = 10, height = 10, append = FALSE)
