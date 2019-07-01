library(hrbrthemes)
library(tidyverse)
library(ggpubr)
library(pwr)

# possibles sample sizes:
Nrange <- seq(2, 400, 2)
ESrange <- seq(0.01, 2, 0.01)
grid <- expand.grid(N = Nrange, ES = ESrange, KEEP.OUT.ATTRS = FALSE)

####################################
# two samples t-test
#############################

power_between <- function(n, es) {
    
    res <- pwr.t.test(
        n = n, d = es, sig.level = 0.05, power = NULL,
        type = "two.sample", alternative = "two.sided"
        )
    
    return(res$power)
    
}

grid %<>%
    # applying the power function each row
    rowwise %>%
    # creating a new column for power
    mutate(power = power_between(N, ES) ) %>%
    ungroup

# identifying the sample size for a power of .90 for the average ES
sample_size <- grid %>%
    #group_by(ES) %>%
    filter(ES == 0.5) %>%
    filter(power >= .90) %>%
    filter(N == min(N) ) %>%
    pull(N)

plotinter <-
    grid %>%
    ggplot(aes(x = ES, y = N) ) +
    geom_raster(aes(fill = power), interpolate = TRUE) +
    scale_fill_gradient(low = "grey10", high = "grey90") +
    geom_line(
        data = . %>% group_by(ES) %>% filter(power >= .90) %>% filter(N == min(N) ),
        colour = "black", size = 1
        ) +
    geom_vline(xintercept = 0.3, lty = 2) +
    geom_vline(xintercept = 0.71, lty = 2) +
    # forcing the origin at zero
    scale_x_continuous(expand = c(0, 0) ) +
    scale_y_continuous(expand = c(0, 0) ) +
    xlab("Effect size") +
    ylab("Sample size") +
    ggtitle("Between-subject") +
    theme_ipsum(
        base_family = "Helvetica",
        base_size = 10, plot_title_size = 12, axis_text_size = 9
        )

############################################
# one sample t-test
################################

power_within <- function(n, es) {
    
    res <- pwr.t.test(
        n = n, d = es, sig.level = 0.05, power = NULL,
        type = "one.sample", alternative = "two.sided"
        )
    
    return(res$power)
    
}

grid <- expand.grid(N = Nrange, ES = ESrange, KEEP.OUT.ATTRS = FALSE) %>%
    # applying the power function each row
    rowwise %>%
    # creating a new column for power
    mutate(power = power_within(N, ES) ) %>%
    ungroup

# identifying the sample size for a power of .90 for the average ES
sample_size <- grid %>%
    filter(ES == 0.5) %>%
    filter(power >= .90) %>%
    filter(N == min(N) ) %>%
    pull(N)

plotintra <-
    grid %>%
    ggplot(aes(x = ES, y = N) ) +
    geom_raster(aes(fill = power), interpolate = TRUE) +
    scale_fill_gradient(low = "grey10", high = "grey90") +
    geom_line(
        data = . %>% group_by(ES) %>% filter(power >= .90) %>% filter(N == min(N) ),
        colour = "black", size = 1
    ) +
    geom_vline(xintercept = 0.3, lty = 2) +
    geom_vline(xintercept = 0.71, lty = 2) +
    # forcing the origin at zero
    scale_x_continuous(expand = c(0, 0) ) +
    scale_y_continuous(expand = c(0, 0) ) +
    xlab("Effect size") +
    ylab("Sample size") +
    ggtitle("Within-subjects") +
    theme_ipsum(
        base_family = "Helvetica",
        base_size = 10, plot_title_size = 12, axis_text_size = 9
        )

ggarrange(plotinter, plotintra, nrow = 1, ncol = 2, legend = "right", common.legend = TRUE)
#ggsave("power.pdf", width = 12, height = 6)
