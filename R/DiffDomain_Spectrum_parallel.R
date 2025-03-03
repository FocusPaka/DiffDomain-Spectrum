#' Parallel codes for TADs detection.
#'
#' @param tadlist_path the TADs need to be detected.
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
#' @param core_num_prop the proportion of core used in parallel calculation.
#' @param seed random seed used in the function of future_map().
#' @param save_path the save path to save the detection information results.
#'
#' @returns all detection results are saved in the save_path file.
#' @export
#'
#' @examples
DiffDomain_Spectrum_parallel <- function(tadlist_path,
                                         fhic0,
                                         fhic1,
                                         min_nbin,
                                         hicnorm,
                                         reso,
                                         prop,
                                         core_num_prop=0.5,
                                         seed=123,
                                         save_path='output.txt'){
    tadlist <- read.table(tadlist_path,header = T)
    future::plan('multisession', workers=round(parallel::detectCores()*core_num_prop))
    y <- furrr::future_map(.x=1:100,     #seq_len(nrow(tadlist)),
                    .f=identifyTADs_HiC,
                    tadlist=tadlist,
                    fhic0 = fhic0,
                    fhic1 = fhic1,
                    min_nbin = min_nbin,
                    hicnorm = hicnorm,
                    reso=reso,prop=prop,
                    .progress = TRUE,
                    .options = furrr::furrr_options(seed = seed))

    combined_df <- do.call(rbind, y)
    row.names(combined_df) <- NULL

    write.table(combined_df, save_path,
                sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
}





