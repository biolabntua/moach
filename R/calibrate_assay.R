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
           Gaussian = mean(x) + qnorm(1 - a) * sd(x),
           Robust = median(x) + qnorm(1 - a) * mad(x),
           Empirical = quantile(x, 1 - a, type = 4))
}

calculateLoD <- function(x, conc, dist = c("Gaussian", "Robust", "Empirical"),
                         a = 0.05, b = 0.05) {
    dist <- match.arg(dist)
    if (b <= 0 || b >= 1) stop("b is a probability, it should be in the [0, 1] interval")
    if (length(x) != length(conc)) stop("x and conc have different lengths")
    LoB <- calculateLoB(x, conc, dist, a)
    low_conc <- (x <= 4 * LoB) & (conc > 0)
    x <- x[low_conc]
    conc <- conc[low_conc]
    if (dist == "Empirical") {
        # assuming constant sd
        Dsb <- aggregate(x ~ conc, FUN = function(y) quantile, probs = b, type = 4)
        Dsb <- mean(Dsb$conc - Dsb$x)
        return(LoB + Dsb)
    }

    df <- length(conc) - length(unique(conc))  # degrees of freedom
    cb <- qnorm(1 - b) / (1 - 1/(4 * df))
    if (dist == "Gaussian") {
        SD <- sqrt(aggregate(x ~ conc, FUN = function(y) length(y) * var(y)) / length(x))
    } else {
        SD <- sqrt(aggregate(x ~ conc, FUN = function(y) length(y) * mad(y) ^ 2) / length(x))
    }
    return(LoB + cb * SD)
}

fit_PLR <- function(dose, response, N = 5, useLog = TRUE, start = NULL, ...) {
    # fit Parametric logistic regression with N parameters, return the fitted object.
    # check: http://www.sciencedirect.com/science/article/pii/S0163725802002231 and
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1472692/
    #
    if (length(dose) != length(response)) stop("dose and response have different lengths.")
    if (N != 4 && N != 5) stop("Only 5 and 4 parametric models are supported")
    df <- data.frame(x = dose, y = response)
    if (useLog) {
        zeros <- df$x == 0
        if (any(zeros)) warning("Zero doses are ignored for log-transformation")
        df <- df[!zeros, , drop = FALSE]
    }

    if (is.null(start)) {
        start <- c(A = min(df$y), B = 1, C = median(df$x), D = max(df$y), E = 1)
    }

    fmla <- switch(as.character(N),
                   "4" = y ~ A + (D - A) / (1 + (x/C) ^ B),
                   "5" = y ~ A + (D - A) / (1 + (x/C) ^ B) ^ E)
    nls(fmla, data = df, start = start)
}

calibration_curve <- function(m) {
    # make a ggplot object ready to plot/save
    NULL
}

calibrate_assay <- function(df) {
    # run calibration pipeline on data frame
    NULL
}
