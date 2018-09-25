#' calcBlockFrame
#'
#' @keywords internal
#'
#' @description  calculates every start and stop position for every block
#'
#' @return a data.frame containing the start and stop positions for each block
calcBlockFrame <-
    function(rowPositions,
             colPositions,
             rowBlockSize,
             colBlockSize) {
        ## construct data frame with all row-and column star-and stop positions
        blocks <- data.frame(
            number = numeric(),
            rowStartPos = numeric(),
            rowStopPos = numeric(),
            colStartPos = numeric(),
            colStopPos = numeric()
        )
        ## counters to choose the correct start - and stop positions within the
        ## loop
        rowCounter <- 1
        colCounter <- 1
        counter <- 1
        ## iteration over rows and columns to calculate the positions
        for (i in 0:(rowPositions[1] - 1)) {
            if (rowCounter < rowPositions[1] - 1) {
                currentRowStart <- i * rowBlockSize + 1
                currentRowStop <- i * rowBlockSize + rowBlockSize
                rowCounter <- rowCounter + 1
            } else {
                if (rowCounter == rowPositions[1] - 1) {
                    currentRowStart <- rowPositions[2]
                    currentRowStop <- rowPositions[4]
                    rowCounter <- rowCounter + 1
                } else {
                    if (rowCounter == rowPositions[1]) {
                        currentRowStart <- rowPositions[3]
                        currentRowStop <- rowPositions[5]
                    }
                }
            }
            for (j in 0:(colPositions[1] - 1)) {
                if (colCounter < colPositions[1] - 1) {
                    currentColStart <- j * colBlockSize + 1
                    currentColStop <- j * colBlockSize + colBlockSize
                    colCounter <- colCounter + 1
                } else {
                    if (colCounter == colPositions[1] - 1) {
                        currentColStart <- colPositions[2]
                        currentColStop <- colPositions[4]
                        colCounter <- colCounter + 1
                    } else {
                        if (colCounter == colPositions[1]) {
                            currentColStart <- colPositions[3]
                            currentColStop <- colPositions[5]
                            colCounter <- 1
                        }
                    }
                }
                blocks <- rbind(
                    blocks,
                    data.frame(
                        number = counter,
                        rowStartPos = currentRowStart,
                        rowStopPos = currentRowStop,
                        colStartPos = currentColStart,
                        colStopPos = currentColStop
                    )
                )
                counter <- counter + 1
            }
        }
        return(blocks)
    }