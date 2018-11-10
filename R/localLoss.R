#' localLoss
#' 
#' @keywords internal
#' 
#' @description Calculates the local loss
#' 
#' @return a list containing two matrices containing the local loss
localLoss2 <- function(L, R, is, js, error_matrix) {
    ## create gradient matrices
    dL <- matrix(0, nrow = nrow(L), ncol = ncol(L))
    dR <- matrix(0, nrow = nrow(R), ncol = ncol(R))
    
    
    ## fill the gradient matrices by calculating the local losses
    for (i in seq_len(length(is))) {
        
        #computation of the local loss of the i-th revealed entry for factorization LR
        x <- error_matrix[is[i], js[i]]
        dL[is[i], ] <- dL[is[i], ] + R[, js[i]] * x
        dR[, js[i]] <- dR[, js[i]] + L[is[i], ] * x
    }
    
    return(list(dL=dL, dR=dR))
}