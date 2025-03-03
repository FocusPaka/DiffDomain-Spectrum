#' The goodness-of-fit test for the reorganized TADs detection.
#'
#' @param x vector to be test.
#' @param y default semicircle.
#' @param N sampled numbers to be used to estimate p-value, default 10000.
#'
#' @returns a list, include the statistic to estimate p-value, the obtained
#' p-value, represent as p.value.
#' @export
#'
#' @examples
#' x <- runif(100, min = -2, max = 2)
#' res <- DiffDomain_Spectrum(x)
#' p_value <- res$p.value
DiffDomain_Spectrum <- function(x, y="semicircle", N=10000){
    DNAME <- deparse(substitute(x))
    x <- x[!is.na(x)]
    n1 <- length(x)
    if (n1 < 8)
        stop("sample size must be greater than 7")

    if (is.character(y)) {
        METHOD <- "One-sample ZA test"
        x <- sort(x)
        distr <- y

        if (distr == "semicircle") {
            p1 <- sapply(x, semicircle_cdf)
            for (i in 1:n1) {
                if (p1[i] == 0) {
                    p1[i] = p1[i] + 1e-10
                }
                if (p1[i] == 1) {
                    p1[i] = p1[i] - 1e-10
                }
            }
        }

        p2 <- 1 - p1
        logp1 <- log(p1)
        logp2 <- log(p2)
        h <- logp1/(n1 - seq(1:n1) + 0.5) + logp2/(seq(1:n1) - 0.5)
        ZA <- -sum(h)
        S <- 0
        if (distr == "semicircle") {
            min_range <- max(c(-2, min(x)))
            max_range <- min(c(2,max(x)))
            if (max_range-min_range<0.01){
                S <- N
            }else{
                for (j in 1:N) {
                    random_numbers <- numeric(0)
                    while (length(random_numbers) < n1) {
                        new_random <- rvs(n1)
                        filtered_random <- new_random[new_random >= min_range &
                                                          new_random <= max_range]
                        random_numbers <- c(random_numbers, filtered_random)
                    }
                    random_numbers <- random_numbers[1:n1]
                    R <- sort(random_numbers)
                    p1 <- sapply(R, semicircle_cdf)
                    for (i in 1:n1) {
                        if (p1[i] == 0) {
                            p1[i] = p1[i] + 1e-10
                        }
                        if (p1[i] == 1) {
                            p1[i] = p1[i] - 1e-10
                        }
                    }
                    p2 <- 1 - p1
                    logp1 <- log(p1)
                    logp2 <- log(p2)
                    h <- logp1/(n1 - seq(1:n1) + 0.5) + logp2/(seq(1:n1) - 0.5)
                    za <- -sum(h)
                    S <- S + (za > ZA)
                }
            }
        }
        p.value <- S/N
    }
    names(ZA) <- "Za"

    RVAL <- list(statistic = ZA, p.value = p.value, alternative = NULL,
                 method = METHOD, data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}
