#' Generate random variables that follow the semi-circle law.
#'
#' @param size The number of generated random variables.
#' @param beta (int, default=1): descriptive integer of the Gaussian ensemble
#' type. For GOE beta=1, for GUE beta=2, for GSE beta=4.
#' @param center (float, default=0.0): center of the distribution.
#' Since the distribution has the shape of a semicircle, the center corresponds
#' to its center.
#' @param sigma (float, default=1.0): scale of the distribution. This value also
#' corresponds to the standard deviation of the random entries of the sampled
#' matrix.
#'
#' @returns The generated random variable.
#' @export
#'
#' @examples
#' x <- rvs(10)
#' x
rvs <- function(size, beta=1, center=0, sigma=1){
    radius <- 2*sqrt(beta)*sigma
    beta_samples <- rbeta(size, shape1 = 1.5, shape2 = 1.5)
    return(center + 2*radius*beta_samples  - radius)
}
