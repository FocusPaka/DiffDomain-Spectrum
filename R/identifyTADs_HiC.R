#' Detect the reorganized TADs between different samples.
#'
#' @param x the index of detected TADs in tadlist.
#' @param tadlist the TADs need to be detected.
#' @param fhic0 the path of .hic file for condition 1.
#' @param fhic1 the path of .hic file for condition 2.
#' @param min_nbin the minimum number of bins for Hi-C contact matrix,
#' default 8.
#' @param hicnorm Normalization to apply. Must be one of NONE/VC/VC_SQRT/KR.
#' VC is vanilla coverage, VC_SQRT is square root of vanilla coverage, and KR
#' is Knight-Ruiz or Balanced normalization.
#' @param reso the resolution of bins.
#' @param prop a threshold value that represent the proportion of zeros in the
#' contact matrix, if the proportion is more than this threshold, then the
#' corresponding rows and columns are removed.
#'
#' @returns a dataframe about the TAD, including the chromosome location, the
#' start and end of TAD, the region of TAD, and the statistics of calculate
#' p-value, the p-value of the TAD, and the size of TAD matrix.
#' @export
#'
#' @examples
#'
identifyTADs_HiC <- function(x,
                            tadlist,
                            fhic0,fhic1,
                            min_nbin,hicnorm,
                            reso,prop){
    chrn <- tadlist[x,][,1]
    start <- tadlist[x,][,2]
    end <- tadlist[x,][,3]
    mat0 <- contact_matrix_from_hic(chrn,start,end,reso,fhic0,hicnorm)
    mat1 <- contact_matrix_from_hic(chrn,start,end,reso,fhic1,hicnorm)

    nbins <- length(makewindow(start,end,reso))
    ind0 <-c() ; ind1 <- c()
    for(i in 1:nrow(mat0)){
        ind0[i] <- sum(is.na(mat0[,i]))
        ind1[i] <- sum(is.na(mat1[,i]))
    }
    ind0 <- ind0 < nbins * prop
    ind1 <- ind1 < nbins * prop
    ind <- ind0 & ind1

    if(sum(ind) >= min_nbin){
        mat0rmna <- mat0[ind,ind]
        mat1rmna <- mat1[ind,ind]

        Diffmat <- mat0rmna / mat1rmna
        Diffmat[which(is.nan(Diffmat))] <- 1
        Diffmatnorm <- normDiffbyMeanSD(Diffmat)

        domname <- paste(sprintf('chr%s',chrn),sprintf(':%s',start),
                         sprintf('-%s',end),sep = '')

        eigv <- eigen(Diffmatnorm/sqrt(nrow(Diffmatnorm)),symmetric = TRUE,
                      only.values = TRUE)$values
        za_stat <- DiffDomain_Spectrum(eigv, N=10000)
        za_p_result <- data.frame(chrn = chrn,start = start,end = end,
                                  region = domname,za=za_stat$statistic,
                                  p = za_stat$p.value, nd=nrow(Diffmatnorm))
        cat(paste0("pvalue_za:",za_p_result$p,"\n"))
    } else {
        domname <- paste(sprintf('chr%s',chrn),sprintf(':%s',start),
                         sprintf('-%s',end),sep = '')
        za_p_result <- data.frame(chrn = chrn,start = start,end = end,
                                  region = domname, za=NA,
                                  p =  NA,nd=NA)
        print('The length of this TAD is too small at this resolution
              to be calculated !')
    }
    return(za_p_result)

}

