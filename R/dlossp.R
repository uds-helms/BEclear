#' dlossp
#' 
#' @keywords internal
#' 
#' @description compute the local loss of the p-th revealed entry for 
#' factorization LR
#' 
#' @return  the result as a list (elements can be accessed with x$Li or x$Ri, 
#' where x is the returned list)
dlossp <- function(L, R, lambda, x, nnzis, nnzjs,i, j, r) {

    
    ## create two length-r vectors dLi and dRj
    ## * dLi[k] contains the gradient of the local loss with
    ## respect to L_ik
    ## * dRi[k] contains the gradient of the local loss with
    ## respect to R_kj
    dLi <- vector(mode="numeric", length=r)
    dRj <- vector(mode="numeric",  length=r)
    
    for(k in seq_len(r)) {
        dLi[k] <- -2 * R[k, j] * (x - L[i, ]%*%R[, j]) + 2 * 
            lambda * L[i, k]/nnzis[i]
        dRj[k] <- -2 * L[i, k] * (x - L[i, ]%*%R[, j]) + 2 * 
            lambda * R[k, j]/nnzjs[j]
    }
    
    return(list(Li=dLi, Rj=dRj))
}