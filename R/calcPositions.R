#' calcPositions
#'
#' @keywords internal
#'
#' @return  returns vector with number of blocks, start - and stop position of
#' last 2 blocks
#'
calcPositions <- function(num, blockSize) {
  if (blockSize <= 0) {
    stop("blockSizes must be greater than 0")
  }

  if (num <= 0) {
    stop("number of entries must be greater than 0")
  }

  rest <- num %% blockSize
  ## normal start and stop position if nrow can be divided by blockSize
  ## without a rest
  startPos1 <- num - (2 * blockSize) - rest + 1
  startPos2 <- num - blockSize - rest + 1
  stopPos1 <- startPos2 - 1
  stopPos2 <- num
  ## special case if matrix has < blockSize rows
  if (num < blockSize) {
    startPos1 <- 0
    startPos2 <- 1
    stopPos1 <- 0
    blocks <- 1
  } else {
    ## num can be divided by blockSize
    if ((num %% blockSize) == 0) {
      blocks <- num / blockSize # done
    } else {
      ## case not dividable by 100
      ## case rest < 75% of block size
      if (rest < floor(blockSize * 0.75)) {
        blocks <- num %/% blockSize # done
      } else {
        ## case rest > 75%
        y <- blockSize + rest
        blocks <- num %/% blockSize + 1 # done
        ## case y dividable by 2
        if (y %% 2 == 0) {
          startPos1 <- num - y + 1
          startPos2 <- startPos1 + (y %/% 2) + 1
          stopPos1 <- startPos2 - 1
        }
        else {
          startPos1 <- num - y + 1
          startPos2 <- startPos1 + (y %/% 2) + 1
          stopPos1 <- startPos2 - 1
        }
      }
    }
  }
  positions <- c(blocks, startPos1, startPos2, stopPos1, stopPos2)
  names(positions) <- c("numBlocks", "startN-1", "startN", "stopN-1", "stopN")
  return(positions)
}
