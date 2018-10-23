#' loss
#' 
#' @keywords internal
#'  
#' @description  computation of the loss of factorization LR
#' 
#' @importFrom Matrix norm
#' 
#' @return The loss calculated
loss <- function(L, R, lambda, is, js, D) {
    lossResult <- 0
    for (i in seq_len(length(is))) {
        lossResult <- lossResult + ((D[is[i], js[i]] -
                                         (L[is[i], ] %*% R[, js[i]]))^2)
    }
    lossResult <- lossResult + lambda * (norm(L, type="F")^2)
    + lambda * (norm(R, type="F")^2)
    
    return(as.numeric(lossResult))
}