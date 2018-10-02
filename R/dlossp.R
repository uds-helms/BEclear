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
    tmpVal <- as.numeric(x - L[i, ] %*% R[, j])
    
    dLi <- -2 * R[, j] * (tmpVal) + 2 * lambda * L[i, ]/nnzis[i]
    dRj <- -2 * L[i,] * (tmpVal) + 2 * lambda * R[, j]/nnzjs[j]
   
    return(list(Li=dLi, Rj=dRj))
}