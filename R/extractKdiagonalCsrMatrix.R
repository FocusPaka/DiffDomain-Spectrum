#' Extract K diagonal data.
#'
#' @param spsCsrMat A symmetric matrix need to be extract K diagonal data.
#'
#' @returns A list, the elements are the values of contact with the same
#' distance in the original symmetric matrix.
#' @export
#'
#' @examples
#' mtx <- matrix(runif(100, min = -1, max = 1), nrow = 10)
#' mtx_Symmetric <- (mtx + t(mtx)) / 2
#' contactBydistance <- extractKdiagonalCsrMatrix(mtx_Symmetric)
#' contactBydistance
extractKdiagonalCsrMatrix <- function(spsCsrMat)
{
    nonZeroIndex <- which(spsCsrMat != 0 | spsCsrMat ==0| is.na(spsCsrMat) |
                              is.infinite(spsCsrMat), arr.ind = TRUE)
    nonZeroIndex_row <- nonZeroIndex[,1]; nonZeroIndex_col <- nonZeroIndex[,2]

    distkey <- c(); distcolIndex <-c(); distrowIndex <-c(); distvalue <-c()

    for(i in 1:length(nonZeroIndex_row)){
        rowIndex <- nonZeroIndex_row[i]
        colIndex <- nonZeroIndex_col[i]
        if(colIndex >= rowIndex){
            distkey <- c(distkey,colIndex - rowIndex)
            distcolIndex <- c(distcolIndex,colIndex)
            distrowIndex <- c(distrowIndex,rowIndex)
            distvalue <- c(distvalue, spsCsrMat[rowIndex,colIndex])
        }
    }
    UpperMatrixIndex <- data.frame(key = distkey,row = distrowIndex,
                                   col = distcolIndex,value = distvalue)
    contactBydistance <- split(UpperMatrixIndex,UpperMatrixIndex$key)
    return(contactBydistance)
}
