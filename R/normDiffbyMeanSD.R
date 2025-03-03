#' Normalize the difference matrix.
#'
#' @param Diffmat A symmetric matrix need to be normalized.
#'
#' @returns A normalized symmetric matrix
#' @export
#'
#' @examples
#' A <- matrix(runif(10^2, min = 0, max = 1), nrow = 10)
#' Diffmat <- (A + t(A)) / 2
#' mtx <- normDiffbyMeanSD(Diffmat)
#' mtx
normDiffbyMeanSD <- function(Diffmat)
{
    options(scipen = 200)
    options(digits = 8)

    D_log <- log(Diffmat)
    contactByDistanceDiff <- extractKdiagonalCsrMatrix(D_log)

    contactMedian <- data.frame(); contactMax <- data.frame();
    contactMean <- data.frame(); contactStd <- data.frame()

    for(k in 1:length(contactByDistanceDiff)){
        val <- contactByDistanceDiff[[k]]$value
        indnan <- is.na(val)
        indinf <- is.infinite((val))

        if(sum(indnan) > 0 & sum(indinf) > 0) {
            ind <- indnan | indinf
            val1 <- val[!ind]
            DataFrameMedian <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                          contactmedian = median(val1))
            contactMedian <- rbind(contactMedian,DataFrameMedian)

            DataFrameMean <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                        contactmean = mean(val1))
            contactMean <- rbind(contactMean,DataFrameMean)

            DataFrameMax <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                       contactmax = max(val1))
            contactMax <- rbind(contactMax,DataFrameMax)

            DataFrameSd <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                      contactstd = sdAdjust(val1))
            contactStd <- rbind(contactStd,DataFrameSd)
        } else if(sum(indnan) > 0 & sum(indinf) == 0) {
            val1 <- val[!indnan]
            if(length(val1) < 1){
                DataFrameMedian <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                              contactmedian = 0)
                contactMedian <- rbind(contactMedian,DataFrameMedian)

                DataFrameMean <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                            contactmean = 0)
                contactMean <- rbind(contactMean,DataFrameMean)
            }else{
                DataFrameMedian <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                              contactmedian = median(val1))
                contactMedian <- rbind(contactMedian,DataFrameMedian)

                DataFrameMean <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                            contactmean = mean(val1))
                contactMean <- rbind(contactMean,DataFrameMean)
            }
            DataFrameSd <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                      contactstd = sdAdjust(val1))
            contactStd <- rbind(contactStd,DataFrameSd)
        } else if (sum(indnan) == 0 & sum(indinf) > 0) {
            val1 <- val[!indinf]
            DataFrameMean <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                        contactmean = mean(val1))
            contactMean <- rbind(contactMean,DataFrameMean)

            DataFrameMax <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                       contactmax = max(val1))
            contactMax <- rbind(contactMax,DataFrameMax)

            DataFrameSd <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                      contactstd = sdAdjust(val1))
            contactStd <- rbind(contactStd,DataFrameSd)
        } else {
            DataFrameMean <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                        contactmean = mean(val))
            contactMean <- rbind(contactMean,DataFrameMean)

            DataFrameSd <- data.frame(key = unique(contactByDistanceDiff[[k]]$key),
                                      contactstd = sdAdjust(val))
            contactStd <- rbind(contactStd,DataFrameSd)
        }
    }

    indnan <- which(is.na(D_log), arr.ind = TRUE)
    if(length(indnan) > 0){
        indr <- indnan[,1]; indc <- indnan[,2]
        for(i in 1:length(indr)){
            k <- abs(indr[i]-indc[i])
            D_log[indr[i],indc[i]] <-
                contactMedian[which(contactMedian$key == k),]$contactmedian
        }
    }

    indinf <- which(is.infinite(D_log), arr.ind = TRUE)
    if(length(indinf) > 0){
        indr <- indinf[,1]; indc <- indinf[,2]
        for(i in 1:length(indr)){
            k <- abs(indr[i]-indc[i])
            D_log[indr[i],indc[i]] <-
                contactMax[which(contactMax$key == k),]$contactmax
        }
    }

    for(i in 1:nrow(contactStd)){
        if(is.na(contactStd$contactstd[i]) | contactStd$contactstd[i] == 0)
            contactStd$contactstd[i] <- 1
    }

    for(i in 1:nrow(D_log)) {
        for(j in 1:nrow(D_log)) {
            k <- abs(i-j)
            D_log[i,j] <-
                (D_log[i,j]-contactMean$contactmean[which(contactMean$key == k)]) /
                (contactStd$contactstd[which(contactStd$key == k)])
        }
    }
    return(D_log)
}
