#' loss
#' 
#' @keywords internal
#'  
#' @description  computation of the loss of factorization LR
#' 
#' @importFrom Matrix norm
#' 
#' @return The loss calculated
loss <- function(L, R, lambda, D) {

    squared_error <- sum((D - (L %*% R))^2)
   
    cost <- lambda * ((norm(L, type="F")^2) + (norm(R, type="F")^2))
    
    
    return(squared_error + cost)
}