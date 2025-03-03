#' The distribution function of the semicircle law.
#'
#' @param x vector of quantiles.
#'
#' @returns A probability value.
#' @export
#'
#' @examples
#' semicircle_cdf(0)
semicircle_cdf <- function(x) {
    stat <- numeric(length(x))
    stat[x>2] <- 1
    stat[x<(-2)] <- 0
    stat[(x>= (-2)) & (x<=2)] <- x[(x>= (-2)) & (x<=2)]*
        sqrt(4-x[(x>= (-2)) & (x<=2)]^2)/(4*pi)+asin(x[(x>= (-2)) & (x<=2)]/2) /
        pi+0.5
    return(stat)
}
