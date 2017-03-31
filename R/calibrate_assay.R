calculateLoB <- function(x, conc = NULL, dist = c("Gaussian", "Robust", "Empirical"),
                         a = 0.05) {
    dist <- match.arg(dist)
    if (a <= 0 || a >= 1) stop("a is a probability measure, it must be in the [0, 1] interval")
    if (!is.null(conc)) {
       if (! length(x) == length(conc)) stop("x and conc have different lengths")
        zeros <- conc == 0
        if (all(!zeros)) stop("No zero concentration was found")
        x <- x[zeros]
    }

    switch(dist,
           Gaussian = qnorm(1 - a, mean(x), sd(x)),
           Robust = qnorm(1 - a, median(x), mad(x)),
           Empirical = quantile(x, 1 - a, type = 4))
}

calculateLoD <- function(x, conc, dist = c("Gaussian", "Robust", "Empirical"),
                         LoB = NULL, b = 0.05, p_var = 0.05, low_conc = NULL) {
    dist <- match.arg(dist)
    if (b <= 0 || b >= 1) stop("b is a probability, it should be in the [0, 1] interval")
    if (length(x) != length(conc)) stop("x and conc have different lengths")
    if (is.null(LoB)) LoB <- calculateLoB(x, conc, dist)
    if (is.null(low_conc)) low_conc <- 4 * LoB

    # prepare data
    AV <- aggregate(x ~ conc, FUN = mean)
    ind_low <- conc %in% AV$conc[(AV$x <= low_conc) & (AV$conc > 0)]
    data <- data.frame(x = x[ind_low], conc = conc[ind_low])
    uconc <- unique(data$conc)  # concentration of samples

    # test for outliers
    var_test <- list(p.value = 1)
    if (length(uconc) == 2) {
        var_test <- with(data, var.test(x[conc == uconc[1]], x[conc == uconc[2]]))
    } else if (length(uconc) > 2) {
        var_test <- outliers::cochran.test(x ~ conc, data)
    }
    if (var_test$p.value < p_var) warning(var_test$alternative)

    # Compute!
    if (dist == "Empirical") {
        # assuming constant sd
        p5 <- aggregate(x ~ conc, data = data, FUN = quantile, probs = b, types = 4)
        AV <- aggregate(x ~ conc, data, FUN = median, na.rm = TRUE)
        Dsb <- mean(AV$x - p5$x)
        # Dsb <- low_conc - quantile(data$x, probs = b, type = 4)
        return(LoB + Dsb)
    }

    norm_test <- aggregate(x ~ conc, data = data, function(y) shapiro.test(y)$p.value)
    if (any(norm_test$x < 0.1)) {
        warning("The normality of the data is disputed by Shapiro & Wilk for some samples")
    }

    df <- nrow(data) - length(uconc)  # degrees of freedom
    cb <- qnorm(1 - b) / (1 - 1/(4 * df))
    if (dist == "Gaussian") {
        SDs <- aggregate(x ~ conc, data = data, FUN = function(y) length(y) * var(y))
        SD <- sqrt(sum(SDs$x) / nrow(data))
    } else {
        SDs <- aggregate(x ~ conc, data = data, FUN = function(y) length(y) * mad(y)^2)
        SD <- sqrt(sum(SDs$x) / nrow(data))
    }
    return(LoB + cb * SD)
}

fit_PLR <- function(dose, response, N = 5, isLog = FALSE, start = NULL, ...) {
    # fit Parametric logistic regression with N parameters, return the fitted object.
    # check: http://www.sciencedirect.com/science/article/pii/S0163725802002231 and
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1472692/
    #
    if (length(dose) != length(response)) stop("dose and response have different lengths.")
    if (N != 4 && N != 5) stop("Only 5 and 4 parametric models are supported")
    df <- data.frame(x = dose, y = response)
    if (!isLog) {
        zeros <- df$x == 0
        if (any(zeros)) warning("Zero doses are ignored for log-transformation")
        df <- df[!zeros, , drop = FALSE]
        df$x <- log(df$x)
    }

    if (is.null(start)) {
        start <- with(df, list(A = min(y), D = max(y), B = 1, E = 1, C = approx(y, x, max(y)/2)$y))
    }

    fmla <- switch(as.character(N),
                   "4" = y ~ D + (A - D) / (1 + exp((x - C) * B)),
                   "5" = y ~ D + (A - D) / (1 + exp((x - C) * B)) ^ E)
    lower <- list(A = 0, D = 0, B = 0, E = 0, C = -Inf)
    nls(fmla, data = df, start = start, lower = lower, nls.control(warnOnly = TRUE))
}

calibration_curve <- function(m) {
    # make a ggplot object ready to plot/save
    NULL
}

calibrate_assay <- function(df) {
    # run calibration pipeline on data frame
    NULL
}
