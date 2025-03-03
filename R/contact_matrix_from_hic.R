#' Extract the contact matrix from a .hic file.
#'
#' @param chrn The chromosome location.
#' @param start The start position of selected contact region.
#' @param end The end position of selected contact region.
#' @param reso The resolution of bins.
#' @param fhic The path of .hic file.
#' @param hicnorm Normalization to apply. Must be one of NONE/VC/VC_SQRT/KR.
#' VC is vanilla coverage, VC_SQRT is square root of vanilla coverage, and KR
#' is Knight-Ruiz or Balanced normalization.
#'
#' @returns The extracted matrix.
#' @export
#'
#' @examples
#'
contact_matrix_from_hic <- function(chrn,start,end,reso,fhic,hicnorm)
{
    options(scipen = 200)
    options(digits = 8)

    domwin <- makewindow(start,end,reso)

    nb <- length(domwin)
    nbs <- 1:nb
    domwin_dict <- data.frame(domwin = domwin,nb = nbs)

    mat <- matrix(data = NA,nrow=nb,ncol=nb)

    region <- stringr::str_c(chrn,domwin[1],domwin[nb],sep = ":")
    print(region)

    if(substring(fhic,stringr::str_length(fhic)-3,stringr::str_length(fhic))  == '.hic'){
        el <- strawr::straw(hicnorm,fhic,region,region,'BP',reso)
        for(i in 1:nrow(el)){
            k <- domwin_dict[which(domwin_dict$domwin == el[i,1]),]$nb
            l <- domwin_dict[which(domwin_dict$domwin == el[i,2]),]$nb
            if(k == l){
                mat[k,l] <- el[i,3]
            }else{
                mat[k,l] <- el[i,3]
                mat[l,k] <- el[i,3]
            }
        }
        return(mat)
    }else{
        sprintf("Sorry, %s dosn't exist" ,fhic)
        mat = NULL
        return(mat)
    }
}
