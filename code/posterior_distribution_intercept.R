library(BEST) # posterior distribution plot

post <- posterior_samples(bmod1)

plotPost(
    post$b_intercept, showCurve = TRUE, showMode = FALSE,
    credMass = NULL, xlab = "", xlim = c(0, 0.6)
    )

intercept <- post$b_intercept
mu <- mean(intercept)
sigma <- sd(intercept)

abline(v = mu + sigma, lty = 2)
abline(v = mu + 2 * sigma, lty = 2)
abline(v = mu + 3 * sigma, lty = 2)
abline(v = mu - sigma, lty = 2)
abline(v = mu - 2 * sigma, lty = 2)
abline(v = mu - 3 * sigma, lty = 2)

library(export) # export to ppt
graph2ppt(file = "/Users/Ladislas/Desktop/plotpost.pptx", width = 10, height = 10, append = TRUE)
