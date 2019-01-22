#' loss
#' @param L a matrix describing the effects of the features
#' @param R a matrix describing the effects of the samples
#' @param D a matrix containing the measured values
#' @param lambda constant that controls the extent of regularization during the
#' gradient descent
#'
#' @description  computation of the loss of factorization LR
#'
#' @importFrom Matrix norm
#'
#' @return a list containing the loss calculated and the error matrix
loss <- function(L, R, lambda, D) {
  error_matrix <- D - (L %*% R)

  squared_error <- sum((error_matrix)^2, na.rm = TRUE)

  cost <- lambda * ((norm(L, type = "F")^2) + (norm(R, type = "F")^2))


  return(list(loss = squared_error + cost, error_matrix = error_matrix))
}
