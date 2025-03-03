#' Calculate the adjust standard deviation.
#'
#' @param numberin A vector to calculate the adjust standard deviation.
#'
#' @returns The adjusted standard deviation.
#' @export
#'
#' @examples
#' x <- rnorm(20)
#' sdAdjust(x)
sdAdjust <- function(numberin)
{
    n <-length(numberin)
    if(n >=2)
        sdA <- sd(numberin)*sqrt(n-1)*sqrt(1/n)
    else
        sdA <- 0
    return((sdA))
}
