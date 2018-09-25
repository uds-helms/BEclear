#' makeBoxplot
#' 
#' @aliases makeBoxplot
#' 
#' @seealso \code{\link{calcScore}}
#' @seealso \code{\link{boxplot}}
#' @seealso \code{\link{correctBatchEffect}}
#' 
#' @keywords hplot
#' 
#' @title produce simple predefined boxplot for methylation data
#' 
#' @description A simple \code{\link{boxplot}} is done with boxes either 
#' separated by batches or by samples and describe the five number summary of 
#' all beta values corresponding to a batch or a sample, respectively. The 
#' batch_ids are shown on the x-axis with a coloring corresponding to the 
#' BEscore.
#' 
#' @details The color code for the batch_ids on the x-axis provides a simple 
#' "traffic light" the user can use to decide if he wants to correct for an 
#' assumed batch effect or not. Green means no batch effect, yellow a possibly 
#' existing not severe batch effect and red stands for an obviously existing 
#' batch effect that should be corrected. The traffic light colors are set 
#' according to the BEscore from the \code{\link{calcScore}} function, values 
#' from 0 to 0.02 are colored in green, from 0.02 to 0.1 in yellow and values 
#' over 0.1 are colored in red.
#' 
#' @param data any matrix filled with beta values, column names have to be 
#' sample_ids corresponding to the ids listed in "samples", row names have to
#' be gene names. 
#' @param samples data frame with two columns, the first column has to contain
#' the sample numbers, the second column has to contain the corresponding batch 
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param score data frame produced by the \code{\link{calcScore}} function.
#' Contains the number of presumably batch affected genes and a BEscore which
#' is needed for the coloring of the batch_ids.
#' @param bySamples should the boxes be separated by samples or not. If not,
#' boxes are separated by the batch_ids.
#' @param col colors for the boxes, refers to the standard \code{\link{boxplot}}
#' R-function. If it is set to "standard", boxes are colored batch-wise (if 
#' separated by samples) or the standard color "yellow" is used (if separated 
#' by batches).
#' @param main main title for the box plot. Default is an empty string.
#' @param xlab label for the x-axis of the box plot. Default is "Batch".
#' @param ylab label for the y-axis of the box plot. Default is "Beta value".
#' @param scoreCol should the batch_ids on the a-axis be colored according to 
#' the BEscore or not? If not, black is used as color for all batch_ids.
#' 
#' @export makeBoxplot
#' @importFrom graphics boxplot mtext par
#' @importFrom methods as
#' @importFrom stats ks.test median rnorm
#' @importFrom utils write.table
#' @usage makeBoxplot(data, samples, score, bySamples=FALSE, col="standard", 
#' main="", xlab="Batch", ylab="Beta value", scoreCol=TRUE)
#' 
#' @return Returns a boxplot on the graphic device with the features explained 
#' above.
#' 
#' @examples 
#' ## Shortly running example. For a more realistic example that takes
#' ## some more time, run the same procedure with the full BEclearData
#' ## dataset.
#' 
#' ## Whole procedure that has to be done to use this function.
#' data(BEclearData)
#' ex.data <- ex.data[31:90,7:26]
#' ex.samples <- ex.samples[7:26,] 
#' 
#' ## Prepare the data for the box plots
#' ## Calculate median difference values and p-values
#' meds <- calcMedians(data=ex.data, samples=ex.samples)
#' pvals <- calcPvalues(data=ex.data, samples=ex.samples) 
#' 
#' ## Summarize p-values and median differences for batch affected genes
#' sum <- calcSummary(medians=meds, pvalues=pvals)
#' 
#' # Calculate the BEscore for the batch_id colorings of the x-axis
#' score <- calcScore(data=ex.data, samples=ex.samples, summary=sum)
#' 
#' ## Simple boxplot for the example data separated by batch
#' makeBoxplot(data=ex.data, samples=ex.samples, score=score, bySamples=FALSE, 
#' main="Some box plot")
#' 
#' ## Simple boxplot for the example data separated by samples
#' makeBoxplot(data=ex.data, samples=ex.samples, score=score, bySamples=TRUE,
#' main="Some box plot")

makeBoxplot <- function(data, samples, score, bySamples=FALSE,
    col="standard", main="", xlab="Batch", ylab="Beta value", scoreCol=TRUE) {
    if (bySamples == FALSE) {
        ## get batch numbers
        batches <- sort(unique(samples$batch_id))
        ## prepare data for the box plot
        boxplotData <-list()
        for (i in seq_len(length(batches)))
        {
            s<-as.character(samples[samples$batch_id == batches[i], 
                "sample_id"])
            boxplotData[[as.character(batches[i])]] <- as.vector(as.matrix(
                data[, s]))
            remove(s)
        }
        annotation <- rep("", length(batches))
        par(las=2)
        if (col == "standard") {
            coloring="yellow"
        }
        else {
            coloring=col
        }
        
        ## make boxplot
        boxplot(boxplotData, col=coloring, main=main, xlab=xlab, ylab=ylab,
            names=annotation, cex.axis=1, xaxt="n")

        if (scoreCol == TRUE) {
            ## get scoring colors for the x-axis labels
            colors <- c()
            for (i in batches) {
                if (score$BEscore[score$batch == i] >= 0.1) {
                    colors <- c(colors, "red")
                } else
                    if (score$BEscore[score$batch == i] >= 0.02 & 
                        score$BEscore[score$batch == i] < 1) {
                            colors <- c(colors, "orange")
                    } else
                        if (score$BEscore[score$batch == i] < 0.02) {
                            colors <- c(colors, "darkgreen")
                        }
            }
        }
        else {
            colors <- rep("black", length(batches))
        }

        ## add colored batch number to x-axis
        for (i in seq_len(length(batches))) {
            mtext(batches[i], side=1, line=1, at=i, las=2, cex=1,
                col=colors[i])
        }
    }
    else {
        ## get batch numbers
        batches <- sort(unique(samples$batch_id))
        ## boxplot separated by sample
        ## sort samples by batch
        samplesNew <- samples
        samplesNew <- samplesNew[order(samplesNew$batch_id), ]
        ## take sample_ids
        barcodes <- samplesNew$sample_id

        ## prepare data
        boxplotData<-list()
        for (i in seq_len(length(barcodes)))
        {
            s <- as.character(samplesNew[samplesNew$sample_id == barcodes[i],
                "sample_id"])
            boxplotData[[as.character(barcodes[i])]] <- as.vector(as.matrix(
                data[, s]))
            remove(s)
        }

        annotation <- rep("", length(barcodes))
        par(las=2)
        ## get colors for the boxes
        batchesForColoring <- sort(samples$batch_id)
        ## store number of samples per batch
        countings <- as.vector(table(batchesForColoring))
        if (col == "standard") {
            coloring <- c()
            colnum <- 1
            for (i in seq_len(length(batches))) {
                if (colnum == 1) {
                    coloring <- c(coloring, rep("yellow", countings[i]))
                    colnum=2
                } else
                    if (colnum == 2) {
                        coloring <- c(coloring, rep("red", countings[i]))
                        colnum=3
                    } else
                        if (colnum == 3) {
                            coloring <- c(coloring, rep("blue", countings[i]))
                            colnum=4
                        } else
                            if (colnum == 4) {
                                coloring <- c(coloring, rep("green",
                                    countings[i]))
                                colnum=1
                            }
            }
        }
        else {
            coloring=col
        }

        ## make the boxplot
        boxplot(boxplotData, main=main, xlab=xlab, ylab=ylab,
            names=annotation, col=coloring, cex.axis=0.8, xaxt="n")

        ## get scoring colors for the x-axis labels
        if (scoreCol == TRUE) {
            colors <- c()
            for (i in batches) {
                if (score$BEscore[score$batch == i] >= 0.1) {
                    colors <- c(colors, "red")
                } else
                    if (score$BEscore[score$batch == i] >= 0.02 & 
                            score$BEscore[score$batch == i] < 1) {
                        colors <- c(colors, "orange")
                    } else
                        if (score$BEscore[score$batch == i] < 0.02) {
                            colors <- c(colors, "darkgreen")
                        }
            }
        }
        else {
            colors <- rep("black", length(batches))
        }

        ## calculate positions to add text for every batch
        positionsText <- c()
        positionCounter <- 0
        for (i in seq_len(length(batches))) {
            if (countings[i] == 1) {
                positionsText <- c(positionsText, positionCounter+1)
            } else
            if (countings[i] == 2) {
                positionsText <- c(positionsText, positionCounter+1)
            } else
            if (countings[i] == 3) {
                positionsText <- c(positionsText, positionCounter+2)
            } else
            ## nearly the mid position of the batch
            if (countings[i] >= 4) {
                positionsText <- c(positionsText, positionCounter+
                    ceiling(countings[i]/2))
            }
        positionCounter <- positionCounter + countings[i]
        }
        
        ## add colored batch number to x-axis
        batchCounter <- 1
        for (i in seq_len(positionCounter)) {
            if (batchCounter == length(batches) + 1) {
                break;
            } else
            ## add number if i is on positionText of batch batchCounter
            if (i == positionsText[batchCounter]) {
                mtext(batches[batchCounter], side=1, line=1, at=i, las=2, 
                    cex=1, col=colors[batchCounter])
                batchCounter <- batchCounter + 1
            } else
            ## add "" if not
            if (i != positionsText[batchCounter]) {
                mtext("", side=1, line=1, at=i, las=2, cex=1)
            }
        }
    }
}
