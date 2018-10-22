#' nnz
#' 
#' @keywords internal
#' 
#' @description counts for every index, the number of non-zero entries
#' 
#' @return a vector n containing the number of non-zero entries per index
nnz <- function(indexes, l) {
    ## array of l entries with 0
    n <- rep(0, l)
    ## get non-zero entries per index
    com <- tapply(indexes, indexes, length)
    ## save number of non-zero entries in n
    n[as.integer(names(com))] <- com
    return(n)
}