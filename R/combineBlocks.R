#' combineBlocks
#' 
#' @keywords internal
#' 
#' @description  load and combine packages
#' 
#' @return the combined matrix consisting of the individual blocks
combineBlocks <- function(blockFrame, rowPos, colPos, dir) {
    D1 <- NULL
    blockCounter <- 1
    for (i in seq_len(rowPos[1])) {
        rowStart <- blockFrame[[2]][[blockCounter]]
        rowStop <- blockFrame[[3]][[blockCounter]]
        for (j in seq_len(colPos[1])) {
            colStart <- blockFrame[[4]][[blockCounter]]
            colStop <- blockFrame[[5]][[blockCounter]]
            ## special case, only one column block
            if (colPos[1] == 1) {
                blockName <- paste("D", blockCounter, sep="")
                filename <- paste(blockName, "row", rowStart, rowStop, 
                                  "col", colStart, colStop, "RData", sep=".")
                file <- paste(dir, filename, sep="/")
                load(file)
                assign(paste(blockName, "row", rowStart, rowStop, "col",
                             colStart, colStop, sep="."), D1)
                frame <- D1
                currentRowGroup <- as.data.frame(frame)
                blockCounter <- blockCounter + 1
            } else {
                if (j == 1) {
                    blockName <- paste("D", blockCounter, sep="")
                    filename <- paste(blockName, "row", rowStart, rowStop,
                                      "col", colStart, colStop, "RData", 
                                      sep=".")
                    file <- paste(dir, filename, sep="/")
                    load(file)
                    assign(paste(blockName, "row", rowStart, rowStop,
                                 "col", colStart, colStop, sep="."), D1)
                    frame <- D1
                    currentRowGroup <- as.data.frame(frame)
                    blockCounter <- blockCounter + 1
                }
                else {
                    blockName <- paste("D", blockCounter, sep="")
                    filename <- paste(blockName, "row", rowStart, rowStop, 
                                      "col", colStart, colStop, "RData", 
                                      sep=".")
                    file <- paste(dir, filename, sep="/")
                    load(file)
                    assign(paste(blockName, "row", rowStart, rowStop,
                                 "col", colStart, colStop, sep="."), D1)
                    frame <- D1
                    currentRowGroup <- cbind(currentRowGroup,
                                             as.data.frame(frame))
                    blockCounter <- blockCounter + 1
                }
            }
        }
        ## special case, only one row block
        if (rowPos[1] == 1) {
            predictedGenes <- currentRowGroup
        } else {
            if (i == 1) {
                predictedGenes <- currentRowGroup
            }
            else {
                predictedGenes <- rbind(predictedGenes, currentRowGroup)
            }
        }
    }
    return(as.matrix(predictedGenes))
}