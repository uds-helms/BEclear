#' dlossp
#' 
#' @keywords internal
#' 
#' @description compute the local loss of the p-th revealed entry for 
#' factorization LR
#' 
#' @return  the result as a list (elements can be accessed with x$Li or x$Ri, 
#' where x is the returned list)
dlossp <- function(Li, Rj, nnzi, nnzj, lambda, error) {

    
    ## create two length-r vectors dLi and dRj
    ## * dLi[k] contains the gradient of the local loss with
    ## respect to L_ik
    ## * dRi[k] contains the gradient of the local loss with
    ## respect to R_kj
    
    #print(Rj * (error))
    
    dLi <- -2 * Rj * (error) + 2 * lambda * Li/nnzi
    dRj <- -2 * Li * (error) + 2 * lambda * Rj/nnzj
   
    return(list(Li=dLi, Rj=dRj))
}