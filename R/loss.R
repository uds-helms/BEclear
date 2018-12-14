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
  error_matrix <- D - (L %*% R)

  squared_error <- sum((error_matrix)^2, na.rm = T)

  cost <- lambda * ((norm(L, type = "F")^2) + (norm(R, type = "F")^2))


  return(list(loss = squared_error + cost, error_matrix = error_matrix))
}
