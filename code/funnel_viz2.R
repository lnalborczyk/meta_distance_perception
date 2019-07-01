viz_funnel2 <- function (x, group = NULL, y_axis = "se", method = "FE", contours = TRUE, 
    sig_contours = TRUE, addev_contours = FALSE, contours_col = "Blues", 
    detail_level = 1, egger = FALSE, trim_and_fill = FALSE, trim_and_fill_side = "left", 
    text_size = 3, point_size = 2, xlab = "Effect", ylab = NULL, 
    group_legend = FALSE, group_legend_title = "", x_trans_function = NULL, 
    x_breaks = NULL) 
{
    if (missing(x)) {
        stop("argument x is missing, with no default.")
    }
    if ("rma" %in% class(x)) {
        es <- as.numeric(x$yi)
        se <- as.numeric(sqrt(x$vi))
        if (method != x$method) {
            message("Note: method argument used differs from input object of class rma.uni (metafor)")
        }
        if (is.null(group) & ncol(x$X) > 1) {
            if (!all(x$X == 1 || x$X == 0) || any(apply(as.matrix(x$X[, 
                -1]), 1, sum) > 1)) {
                stop("Can not deal with metafor output object with continuous and/or more than one categorical moderator variable(s).")
            }
            no.levels <- ncol(x$X) - 1
            group <- factor(apply(as.matrix(x$X[, -1]) * rep(1:no.levels, 
                each = length(es)), 1, sum))
        }
    }
    else {
        if ((is.data.frame(x) || is.matrix(x)) && ncol(x) >= 
                2) {
            if (sum(is.na(x[, 1])) != 0 || sum(is.na(x[, 2])) != 
                    0) {
                warning("The effect sizes or standard errors contain missing values, only complete cases are used.")
                if (!is.null(group)) {
                    group <- group[stats::complete.cases(x)]
                }
                x <- x[stats::complete.cases(x), ]
            }
            if (!is.numeric(x[, 1]) || !is.numeric(x[, 2])) {
                stop("Input argument has to be numeric; see help(viz_funnel) for details.")
            }
            if (!all(x[, 2] > 0)) {
                stop("Non-positive standard errors supplied")
            }
            es <- x[, 1]
            se <- x[, 2]
        }
        else {
            stop("Unknown input argument; see help(viz_funnel) for details.")
        }
    }
    if (!is.null(group) && !is.factor(group)) {
        group <- as.factor(group)
    }
    if (!is.null(group) && (length(group) != length(es))) {
        warning("length of supplied group vector does not correspond to the number of studies; group argument is ignored")
        group <- NULL
    }
    k <- length(es)
    summary_es <- metafor::rma.uni(yi = es, sei = se, method = method)$b[[1]]
    summary_se <- sqrt(metafor::rma.uni(yi = es, sei = se, method = method)$vb[[1]])
    summary_tau2 <- metafor::rma.uni(yi = es, sei = se, method = method)$tau2
    if (is.null(group)) {
        plotdata <- data.frame(es, se)
    }
    else {
        plotdata <- data.frame(es, se, group)
    }
    if (!(contours_col %in% c("Blues", "Greys", "Oranges", "Greens", 
        "Reds", "Purples"))) {
        warning("Supported arguments for contours_col are Blues, Greys, Oranges, Greens, Reds, and Purples. Blues is used.")
        contours_col <- "Blues"
    }
    col <- RColorBrewer::brewer.pal(n = 9, name = contours_col)
    if (detail_level < 0.1) {
        detail_level <- 0.1
        warning("Argument detail_level too low. Set to minimum value (0.1)")
    }
    if (detail_level > 10) {
        detail_level <- 10
        warning("Argument detail_level too high. Set to minimum value (10)")
    }
    min_x <- min(plotdata$es)
    max_x <- max(plotdata$es)
    if (trim_and_fill == TRUE) {
        trimnfill <- function(es, se, group = NULL, side = "left") {
            if (side == "right") {
                es <- -es
            }
            if (side != "right" && side != "left") {
                stop("trim_and_fill_side argument must be either left or right")
            }
            mean_func <- function(es, se) {
                metafor::rma.uni(yi = es, sei = se, method = method)$b[1]
            }
            k0_func <- function(es, se, summary_es) {
                n <- length(es)
                Tn <- sum(rank(abs(es - summary_es))[sign(es - 
                        summary_es) > 0])
                round(max((4 * Tn - n * (n + 1))/(2 * n - 1), 
                    0), 0)
            }
            summary_es_init <- mean_func(es, se)
            k0 <- k0_func(es = es, se = se, summary_es = summary_es_init)
            eps <- 1
            iter <- 0
            while (eps > 0.01 || iter < 20) {
                iter <- iter + 1
                es_ord <- es[order(es, decreasing = T)]
                se_ord <- se[order(es, decreasing = T)]
                if (k0 > 0) {
                    es_ord <- es_ord[-(1:k0)]
                    se_ord <- se_ord[-(1:k0)]
                }
                summary_es_new <- mean_func(es_ord, se_ord)
                k0 <- k0_func(es = es, se = se, summary_es = summary_es_new)
                eps <- abs(summary_es_init - summary_es_new)
                summary_es_init <- summary_es_new
            }
            if (iter == 19) {
                warning("Trim and fill algorithm did not converge after 10 iterations")
            }
            if (k0 > 0) {
                es_ord <- es[order(es, decreasing = T)]
                se_ord <- se[order(es, decreasing = T)]
                if (!is.null(group)) {
                    group_ord <- group[order(es, decreasing = T)]
                    group_fill <- group_ord[1:k0]
                }
                if (side == "right") {
                    es_fill <- -(summary_es_new + (summary_es_new - 
                            es_ord[1:k0]))
                    summary_es_init <- -summary_es_init
                }
                else {
                    es_fill <- summary_es_new + (summary_es_new - 
                            es_ord[1:k0])
                }
                se_fill <- se_ord[1:k0]
                if (is.null(group)) {
                    data.frame(es_fill, se_fill, summary_es_init)
                }
                else {
                    data.frame(es_fill, se_fill, group_fill, summary_es_init)
                }
            }
            else {
                if (is.null(group)) {
                    data.frame(es_fill = NULL, se_fill = NULL, 
                        summary_es_init = NULL)
                }
                else {
                    data.frame(es_fill = NULL, se_fill = NULL, 
                        group_fill = NULL, summary_es_init = NULL)
                }
            }
        }
        side <- trim_and_fill_side
        if (is.null(group)) {
            tnfdata <- trimnfill(es, se, side = side)
        }
        else {
            tnfdata <- trimnfill(es, se, group, side = side)
        }
        if (nrow(tnfdata) > 0) {
            if (is.null(group)) {
                names(tnfdata) <- c("es", "se", "tnf_summary")
            }
            else {
                names(tnfdata) <- c("es", "se", "group", "tnf_summary")
            }
            min_x <- min(c(min_x, min(tnfdata$es)))
            max_x <- max(c(max_x, max(tnfdata$es)))
        }
        else {
            trim_and_fill <- FALSE
        }
    }
    if (method == "DL" && addev_contours == TRUE) {
        rem_dl <- function(es, se) {
            summary_es_FEM <- sum((1/se^2) * es)/sum(1/se^2)
            n <- length(es)
            if (n == 1) {
                t2 <- 0
            }
            else {
                Q <- sum((1/se^2) * (es - summary_es_FEM)^2)
                t2 <- max(c(0, (Q - (n - 1))/(sum(1/se^2) - sum((1/se^2)^2)/sum(1/se^2))))
            }
            w <- 1/(se^2 + t2)
            c(sum(w * es)/sum(w), sqrt(1/sum(w)))
        }
    }
    if (y_axis == "se") {
        plotdata$y <- se
        max_se <- max(se) + ifelse(length(se) > 1, diff(range(se)) * 
                0.1, max(se) * 0.1)
        y_limit <- c(0, max_se)
        if (is.null(ylab)) {
            ylab <- "Standard Error"
        }
        if (trim_and_fill == TRUE) {
            tnfdata$y <- tnfdata$se
        }
        if (sig_contours == TRUE) {
            sig_funneldata <- data.frame(x = c(-stats::qnorm(0.975) * 
                    max_se, 0, stats::qnorm(0.975) * max_se, stats::qnorm(0.995) * 
                    max_se, 0, -stats::qnorm(0.995) * max_se), y = c(max_se, 
                        0, max_se, max_se, 0, max_se))
            min_x <- min(c(min_x, min(sig_funneldata$x)))
            max_x <- max(c(max_x, max(sig_funneldata$x)))
        }
        if (contours == TRUE) {
            funneldata <- data.frame(
                x = c(
                    summary_es - stats::qnorm(0.975) * sqrt(max_se^2 + summary_tau2),
                    summary_es,# - 1.96 * sqrt(summary_tau2),
                    summary_es,# + 1.96 * sqrt(summary_tau2),
                    summary_es + stats::qnorm(0.975) * sqrt(max_se^2 + summary_tau2)
                    ),
                
                y = c(max_se, 0, 0, max_se)
                
                )
            
            min_x <- min(c(min_x, min(funneldata$x)))
            max_x <- max(c(max_x, max(funneldata$x)))
        }
        if (egger == TRUE) {
            plotdata <- data.frame(plotdata, z = (plotdata$es)/plotdata$se)
            plotdata <- data.frame(plotdata, prec = 1/plotdata$se)
            radial_intercept <- stats::coef(stats::lm(z ~ prec, 
                data = plotdata))[1]
            radial_slope <- stats::coef(stats::lm(z ~ prec, data = plotdata))[2]
            eggerdata <- data.frame(intercept = radial_slope/radial_intercept, 
                slope = -1/radial_intercept)
        }
    }
    else {
        if (y_axis == "precision") {
            plotdata$y <- 1/se
            max_y <- max(1/se) + ifelse(length(se) > 1, diff(range(1/se)) * 
                    0.05, 1/se * 0.05)
            min_y <- min(1/se) - ifelse(length(se) > 1, diff(range(1/se)) * 
                    0.05, 1/se * 0.05)
            if (is.null(ylab)) {
                ylab <- "Precision (1/SE)"
            }
            if (trim_and_fill == TRUE) {
                tnfdata$y <- 1/tnfdata$se
            }
            if (sig_contours == TRUE) {
                n_support <- 200 * detail_level
                prec <- seq(from = min_y, to = max_y, length.out = n_support)
                x_prec_0.05 <- stats::qnorm(0.975) * (1/prec)
                x_prec_0.01 <- stats::qnorm(0.995) * (1/prec)
                sig_funneldata <- data.frame(x = c(-x_prec_0.01, 
                    rev(x_prec_0.01), x_prec_0.05, rev(-x_prec_0.05)), 
                    y = c(prec, rev(prec), prec, rev(prec)))
                min_x <- min(c(min_x, min(sig_funneldata$x)))
                max_x <- max(c(max_x, max(sig_funneldata$x)))
            }
            if (contours == TRUE) {
                n_support <- 200 * detail_level
                prec <- seq(from = min_y, to = max_y, length.out = n_support)
                x_prec <- stats::qnorm(0.975) * sqrt((1/prec)^2 + 
                        summary_tau2)
                funneldata <- data.frame(x = rep(summary_es, 
                    times = n_support * 2) + c(-x_prec, rev(x_prec)), 
                    y = c(prec, rev(prec)))
                min_x <- min(c(min_x, min(funneldata$x)))
                max_x <- max(c(max_x, max(funneldata$x)))
            }
            if (egger == TRUE) {
                warning("Note: egger = TRUE ignored: Egger's regression line can only be plotted for y_axis = se")
            }
            y_limit <- c(min_y, max_y)
        }
        else {
            stop("y_axis argument must be either se or precision")
        }
    }
    x_limit <- c(min_x - diff(c(min_x, max_x)) * 0.05, max_x + 
            diff(c(min_x, max_x)) * 0.05)
    if (addev_contours == TRUE) {
        if (y_axis == "se") {
            y_range <- c(0.001, max_se + diff(range(y_limit)) * 
                    0.2)
            x_range <- c(min_x - diff(range(x_limit)) * 0.2, 
                max_x + diff(range(x_limit)) * 0.2)
            step <- abs(summary_es - x_range[1])/(150 * detail_level - 
                    1)
            x_add <- c(seq(from = x_range[1], to = summary_es, 
                length.out = 150 * detail_level), seq(from = summary_es + 
                        step, to = x_range[2], by = step))
            y_add <- seq(from = y_range[1], to = y_range[2], 
                length.out = length(x_add))
        }
        else {
            y_range <- c(max_y + diff(range(x_limit)) * 0.2, 
                min_y - diff(range(x_limit)) * 0.2)
            x_range <- c(min_x - diff(range(x_limit)) * 0.2, 
                max_x + diff(range(x_limit)) * 0.2)
            step <- abs(summary_es - x_range[1])/(150 * detail_level - 
                    1)
            x_add <- c(seq(from = x_range[1], to = summary_es, 
                length.out = 150 * detail_level), seq(from = summary_es + 
                        step, to = x_range[2], by = step))
            y_add <- 1/seq(from = y_range[1], to = y_range[2], 
                length.out = length(x_add))
        }
        study_grid <- expand.grid(x_add, y_add)
        names(study_grid) <- c("x_add", "y_add")
        addev_data <- apply(study_grid, 1, function(x) {
            if (method == "FE") {
                M_new <- sum((1/c(se, x[2])^2) * c(es, x[1]))/sum(1/c(se, 
                    x[2])^2)
                Mse_new <- sqrt(1/sum(1/c(se, x[2])^2))
                p.val <- stats::pnorm(M_new/Mse_new)
                c(M_new, p.val)
            }
            else {
                if (method == "DL") {
                    res_dl <- rem_dl(es = c(es, x[1]), se = c(se, 
                        x[2]))
                    M_new <- res_dl[1]
                    p.val <- stats::pnorm(res_dl[1]/res_dl[2])
                    c(M_new, p.val)
                }
                else {
                    mod <- metafor::rma.uni(yi = c(es, x[1]), sei = c(se, 
                        x[2]), method = method, control = list(stepadj = 0.5, 
                            maxiter = 1000))
                    p.val <- stats::pnorm(mod$z)
                    M_new <- mod$b[[1]]
                    c(M_new, p.val)
                }
            }
        })
        addev_data <- t(addev_data)
        addev_data <- data.frame(study_grid, M = addev_data[, 
            1], sig_group = factor(ifelse(addev_data[, 2] < 0.025, 
                "sig.neg. ", ifelse(addev_data[, 2] > 0.975, "sig.pos. ", 
                    "not sig. ")), levels = c("sig.neg. ", "not sig. ", 
                        "sig.pos. ")))
        addev_data <- addev_data[order(addev_data$x_add, decreasing = F), 
            ]
        if (y_axis == "precision") {
            addev_data$y_add <- 1/addev_data$y_add
        }
    }
    if (!is.null(x_trans_function) && !is.function(x_trans_function)) {
        warning("Argument x_trans_function must be a function; input ignored.")
        x_trans_function <- NULL
    }
    y <- NULL
    sig_group <- NULL
    x.01 <- NULL
    x.05 <- NULL
    tnf_summary <- NULL
    intercept <- NULL
    slope <- NULL
    p <- ggplot(data = plotdata, aes(x = es, y = y))
    if (addev_contours == TRUE) {
        p <- p + geom_raster(data = addev_data, aes(x = x_add, 
            y = y_add, fill = sig_group), alpha = 0.4) + scale_fill_manual(name = "", 
                values = c(col[9], col[1], col[4]), drop = FALSE)
    }
    if (sig_contours == TRUE && y_axis == "se") {
        p <- p + geom_polygon(data = sig_funneldata, aes(x = x, 
            y = y), fill = col[9], alpha = 0.6) + geom_path(data = sig_funneldata, 
                aes(x = x, y = y))
    }
    else {
        if (sig_contours == TRUE && y_axis == "precision") {
            p <- p + geom_polygon(data = sig_funneldata, aes(x = x, 
                y = y), fill = col[9], alpha = 0.6) + geom_path(data = sig_funneldata, 
                    aes(x = x, y = y))
        }
    }
    if (contours == TRUE) {
        p <- p + geom_path(data = funneldata, aes(x = x, y = y)) + 
            geom_vline(xintercept = summary_es, linetype = 2)
    }
    if (y_axis == "se") {
        p <- p + scale_y_reverse(name = ylab)
    }
    else {
        if (y_axis == "precision") {
            p <- p + scale_y_continuous(name = ylab)
        }
    }
    if (trim_and_fill == TRUE) {
        if (dim(tnfdata)[1] > 0) {
            if (is.null(group)) {
                p <- p + geom_point(data = tnfdata, aes(x = es, 
                    y = y), size = point_size, col = "black", alpha = 1)
            }
            else {
                p <- p + geom_point(data = tnfdata, aes(x = es, 
                    y = y, shape = group), size = point_size, col = "black", 
                    alpha = 1)
            }
            if (contours == TRUE) {
                p <- p + geom_vline(data = tnfdata, aes(xintercept = tnf_summary), 
                    lty = "dashed")
            }
        }
    }
    if (is.null(group)) {
        p <- p + geom_point(size = point_size, fill = "white", 
            shape = 21, col = "black", alpha = 1)
    }
    else {
        p <- p + geom_point(aes(col = group, shape = group), 
            size = point_size, alpha = 1)
    }
    if (egger == TRUE && y_axis == "se") {
        p <- p + geom_abline(data = eggerdata, aes(intercept = intercept, 
            slope = slope), lty = "dashed", lwd = 1, color = "firebrick")
    }
    if (!is.null(x_trans_function)) {
        if (is.null(x_breaks)) {
            p <- p + scale_x_continuous(name = xlab, labels = function(x) {
                round(x_trans_function(x), 3)
            })
        }
        else {
            p <- p + scale_x_continuous(name = xlab, labels = function(x) {
                round(x_trans_function(x), 3)
            }, breaks = x_breaks)
        }
    }
    else {
        if (is.null(x_breaks)) {
            p <- p + scale_x_continuous(name = xlab)
        }
        else {
            p <- p + scale_x_continuous(breaks = x_breaks, name = xlab)
        }
    }
    p <- p + coord_cartesian(xlim = x_limit, ylim = y_limit, 
        expand = F) + scale_shape_manual(values = 15:19, name = group_legend_title) + 
        scale_color_brewer(name = group_legend_title, palette = "Set1", 
            type = "qual")
    if (group_legend == FALSE) {
        p <- p + guides(color = "none", shape = "none")
    }
    if (addev_contours == TRUE) {
        legend.key <- element_rect(color = "black")
    }
    else {
        legend.key <- element_rect(color = "white")
    }
    p <- p + theme_bw() + theme(text = element_text(size = 1/0.352777778 * 
            text_size), legend.position = "bottom", legend.key = legend.key, 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    p
}
