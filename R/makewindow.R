#' Calculate the positions of window in a region.
#'
#' @param start The start position of a region.
#' @param end The end position of a region.
#' @param reso The resolution of bin.
#'
#' @returns The positions of window.
#' @export
#'
#' @examples
#' wins <- makewindow(1, 100000, 10000)
#' wins
makewindow <- function(start,end,reso)
{
    k1 <- floor(start / reso)
    k2 <- ceiling(end / reso)
    wins <- seq(from = k1*reso, to = k2*reso, by = reso)
    return(wins)
}

