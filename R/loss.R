#' loss
#' 
#' @keywords internal
#'  
#' @description  computation of the loss of factorization LR
#' 
#' @return The loss calculated
loss <- function(L, R, lambda, N, is, js, D) {
  lossResult <- 0
  for (i in seq_len(N)) {
    lossResult <- lossResult + ((D[is[i], js[i]] -
                                   (L[is[i], ]%*%R[, js[i]]))^2)
  }
  lossResult <- lossResult + lambda * (norm(L, type="F")^2)
  + lambda * (norm(R, type="F")^2)
  
  return(as.numeric(lossResult))
}