#' imputeMissingDataForBlock
#'
#' @import futile.logger
#' @import Matrix
#' @importFrom stats rnorm
#'
#' @param fixedSeed determines if they seed should be fixed, which is important
#' for testing
#'
#' @keywords internal
#'
#' @return number of the block processed
imputeMissingDataForBlock <- function(block, dir, epochs, lambda = 1,
                                      gamma = 0.01, r = 10, fixedSeed = TRUE) {
  blockNr <- block$blockNr
  D <- block$block
  flog.info(paste("Impute missing data for block", blockNr, "of", block$total))


  ## check if NA values are contained in the block
  if (any(is.na(D))) {


    ## set D as matrix
    dat <- as.matrix(D)

    ## set NA data to 0
    D <- dat
    D[is.na(D)] <- 0
    D <- as(as.matrix(D), "dgCMatrix")

    Dsummary <- summary(D)
    m <- nrow(D) # number of rows
    n <- ncol(D) # number of columns
    is <- Dsummary$i # row of each revealed entry
    js <- Dsummary$j # column of each revealed entry
    xs <- Dsummary$x # value of each revealed entry
    N <- length(is) # number of revealed entries


    if (fixedSeed) {
      set.seed(1, kind = "Mersenne-Twister")
    }

    L <- matrix(rnorm(m * r), m, r) / sqrt(r)
    R <- matrix(rnorm(r * n), r, n) / sqrt(r)




    ## run LFM
    resultGd <- runGradientDescent(
      L = L, R = R, lambda = lambda,
      epochs = epochs, gamma = gamma,
      blockNr = blockNr, is = is, js = js,
      D = dat, r = r
    )

    dataTemp <- block$block
    ## Predicted matrix

    D1 <- resultGd$L %*% resultGd$R
    ## Replace predicted values with known, keep predicted unknown
    ## values
    for (i1 in seq_len(nrow(dataTemp)))
      for (j1 in seq_len(ncol(dataTemp)))
      {
        if (!is.na(dataTemp[i1, j1])) {
          D1[i1, j1] <- dataTemp[i1, j1]
        }
      }
  }

  ## no NA values contained in the block - keep original values
  else {
    flog.debug(paste(
      "Block", blockNr,
      "has no missing values. Original values are kept"
    ))
    D1 <- D
  }

  D1 <- as.matrix(D1)
  colnames(D1) <- colnames(D)
  rownames(D1) <- rownames(D)

  filename <- paste("D", blockNr, ".RData", sep = "")
  save(D1, file = paste(dir, filename, sep = "/"))

  
  return(list(blockNr = block$blockNr, block = D1, total = block$total))
}
