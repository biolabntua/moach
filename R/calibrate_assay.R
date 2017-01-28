calculateLoB <- function(x, conc = NULL, dist = c("Gaussian", "Robust", "Empirical"), a = 0.05) {
    dist <- match.arg(dist)
    if (!is.null(conc)) x <- x[conc == 0]
    # TODO: add option for student dist?
    if (dist == "Gaussian") {
        return(mean(x) + qnorm(1 - a) * sd(x))
    } else if (dist == "Robust") {  # estimate mu, sigma with median and mad
        # check here: http://stats.stackexchange.com/questions/67337/mad-in-relation-to-95-confidence
        madx <- mad(x, constant = 1)
        x <- x[(x - madx > 4) & (x + madx < 4)]  # trim outliers
        return(median(x) + qnorm(1 - a) * mad(x))
    } else if (dist == "Empirical") {
        return(quantile(x, 1 - a))
    }
}

calculateLoD <- function(x, conc, dist = c("Gaussian", "Robust", "Empirical"), a = 0.05) {
    dist <- match.arg(dist)
}

fit_PLR <- function(x, y, N = 4, useLog = TRUE) {
    # fit Parametric logistic regression with N parameters, return the fitted object.
    NULL
}

calibration_curve <- function(m) {
    # make a ggplot object ready to plot/save
    NULL
}

calibrate_assay <- function(df) {
    # run calibration pipeline on data frame
    NULL
}
