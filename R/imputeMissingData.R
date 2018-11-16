#'imputeMissingData 
#'
#'@aliases imputeMissingData 
#'
#'@title Matrix prediction method using a Latent Factor Model
#'
#'@description This function predicts the missing entries of an input matrix 
#'(NA values) through the use of a Latent Factor Model. You can run the 
#'function also in parallel mode and split up the matrix to a certain number of 
#'smaller matrices to speed up the prediction process. If you set the 
#'rowBlockSize and colBlockSize both to 0, the function is running on the whole
#'matrix. Take a look at the details section for some deeper information about 
#'this. The default parameters are chosen with the intention to make an 
#'accurate prediction within an affordable time interval.
#'
#'@details The method used to predict the missing entries in the matrix is 
#'called "latent factor model". In the following sections, the method itself is 
#'described as well as the correct usage of the parameters. The parameters are 
#'described in the same order as they appear in the usage section.\cr
#'The method originally stems from recommender systems where the goal is to 
#'predict user ratings of products. It is based on matrix factorization and uses
#'a discrete gradient descent (GDE) algorithm that stepwise predicts two 
#'matrices L and R with matching dimensions to the input matrix. These two 
#'matrices are initialized with random numbers and stepwise adjusted towards 
#'the values of the input matrix through the GDE algorithm. After every 
#'adjustment step, the global loss is calculated and the parameters used for 
#'the adjustment are possibly also adjusted so that the global loss is getting 
#'minimized and the prediction is getting accurate. After a predefined number 
#'of steps (called epochs) are executed by the GDE algorithm, the predicted 
#'matrix is calculated by matrix multiplication of L and R. Finally, all 
#'missing values in the input matrix are replaced with the values from the 
#'predicted matrix and the already known values from the input matrix are 
#'maintained. The completed input matrix is then returned at the end.\cr
#' Description of the parameters:
#' \itemize{
#' \item data: simply the input matrix with missing values set to NA
#' \item rowBlockSize and colBlockSize: Here you can define the dimensions of 
#' the smaller matrices, the input matrix is divided into if the function is 
#' working in parallel mode. For details about these so called blocks, see the 
#' section "About the blocks" below.
#' \item epochs: Defines the number of steps the gradient descent algorithm 
#' performs until the prediction ends. Note that the higher this number is, the 
#' more precisely is the prediction and the more time is needed to perform the 
#' prediction. If the step size is too small, the prediction would not be very 
#' good. We suggest to use a step size of 50 since we did not get better 
#' predictions if we took higher step sizes during our testing process.
#' }
#' About the blocks:
#' You have the possibility to change the size of the blocks in which the input 
#' matrix can be divided. if you choose e.g. the rowBlockSize = 50 and the 
#' colBlockSize = 60 your matrix will be cut into smaller matrices of the size 
#' approximately 50x60. Note that this splitting algorithm works with every 
#' possible matrix size! If both size parameters do not fit to the dimensions of
#' the input matrix, the remaining rows and columns of the input matrix are 
#' distributed over some blocks, so that the block sizes are roughly of the same
#' size. All blocks are saved at the specified directory after the processing 
#' of a block has been done within an RData file. These RData files are 
#' continuously numbered and contain the row and column start and stop 
#' positions in their name. Next, these blocks are assembled into the returned 
#' matrix and this matrix is saved in the specified directory. Finally, single 
#' blocks are deleted. To see how this is done, simply run the example at the 
#' end of this documentation. We suggest to use the block size of 60 (default) 
#' but you can also use any other block size, as far as it is bigger than the 
#' number of samples in the biggest batch. This avoids having an entire row of 
#' NA values in a block which leads to a crash of the imputeMissingData method.
#' In order  to process the complete matrix without dividing into blocks, 
#' specify rowBlockSize = 0 and colBlockSize = 0. But if the input matrix is 
#' large (more than 200x200), it is not recommended due to exponential increase
#' of computation time required.\cr
#' Note that the size of the blocks affect the prediction accuracy. In case of 
#' very small blocks, the information obtained from neighbor entries is not 
#' sufficient. Thus, the larger the size of the block is, the more accurately 
#' those entries are predicted. Default size 60 is enough to have accurate 
#' prediction in a reasonable amount of time.
#' 
#' @return Returns a data matrix with the same dimensions as well as same row 
#' and column names as the input matrix. According to the "outputFormat"
#' parameter, either a .RData file containing only the returned matrix or a 
#' tab-delimited .txt file containing the content of the returned matrix is 
#' saved in the specified directory.
#' 
#' @references \insertRef{Akulenko2016}{BEclear}
#' @references \insertRef{Koren2009}{BEclear}
#' @references \insertRef{Candes2009}{BEclear}
#' 
#' @param data any matrix filled e.g. with beta values. The missing entries you
#' want to predict have to be set to NA
#' @param rowBlockSize the number of rows that is used in a block if the 
#' function is run in parallel mode and/or not on the whole matrix. Set this and
#' the "colBlockSize" parameter to 0 if you want to run the function on the 
#' whole input matrix. We suggest to use a block size of 60 but you can also use
#' any other block size, but the size has to be bigger than the number of 
#' samples in the biggest batch. Look at the details section for more 
#' information about this feature.
#' @param colBlockSize the number of columns that is used in a block if the 
#' function is run in parallel mode and/or not on the whole matrix. Set this, 
#' and the "rowBlockSize" parameter to 0 if you want to run the function on the 
#' whole input matrix. We suggest to use a block size of 60 but you can also use
#' any other block size, but the size has to be bigger than the number of 
#' samples in the biggest batch. Look at the details section for more 
#' information about this feature.
#' @param epochs the number of iterations used in the gradient descent algorithm
#' to predict the missing entries in the data matrix.
#' @param lambda constant that controls the extent of regularization during the 
#' gradient descent
#' @param gamma constant that controls the extent of the shift of parameters 
#' during the gradient descent
#' @param r length of the second dimension of variable matrices R and L
#' @param outputFormat you can choose if the finally returned data matrix should
#' be saved as an .RData file or as a tab-delimited .txt file in the specified 
#' directory. Allowed values are "RData" and "txt".
#' @param dir set the path to a directory the predicted matrix should be 
#' stored. The current working directory is defined as default parameter.
#' @param BPPARAM An instance of the 
#' \code{\link[BiocParallel]{BiocParallelParam-class}} that determines how to 
#' parallelisation of the functions will be evaluated.
#' @param matrixOfOnes instead of starting with a random matrix, start from a matrix
#' of ones. For testing purposes only!
#' 
#' @export imputeMissingData 
#' @import BiocParallel
#' @import futile.logger
#' @usage imputeMissingData(data, rowBlockSize=60,  colBlockSize=60, epochs=50, 
#' lambda = 1, gamma = 0.01, r = 10, outputFormat="RData", dir=getwd(), 
#' BPPARAM=bpparam(), matrixOfOnes = FALSE)
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
#' # Calculate median difference values and p-values
#' library(data.table)
#' samples <- data.table(ex.samples)
#' data <- data.table(feature=rownames(ex.data), ex.data)
#' data <- melt(data = data, id.vars = "feature", variable.name = "sample", 
#'             value.name = "beta.value")
#' setkey(data, "feature", "sample")
#' meds <- calcMedians(data=data, samples=samples)
#' pvals <- calcPvalues(data=data, samples=samples)
#' # Summarize p-values and median differences for batch affected genes
#' sum <- calcSummary(medians=meds, pvalues=pvals)
#' 
#' # Set entries defined by the summary to NA
#' clearedMatrix <- clearBEgenes(data=ex.data, samples=ex.samples, summary=sum)
#' 
#' # Predict the missing entries with standard values, row- and block sizes are
#' # just set to 10 to get a short runtime. To use these parameters, either use
#' # the default values or please note the description in the details section
#' # above
#' predicted <- imputeMissingData (data=clearedMatrix, rowBlockSize=10, 
#' colBlockSize=10)

imputeMissingData <- function(data, rowBlockSize=60, colBlockSize=60, epochs=50,
                              lambda = 1, gamma = 0.01, r = 10,
                              outputFormat="RData", dir=getwd(), 
                              BPPARAM=bpparam(), matrixOfOnes = FALSE) {
    
    flog.info("Starting the imputation of missing values.")
    flog.info("This might take a while.")
    D1 <- NULL
    if (epochs<= 0) {
        stop('number of epochs has to be greater than 0')
    } 
    
    ## run BEclear 
    flog.info("BEclear is startet in parallel mode:")
    flog.info(paste("block size:", rowBlockSize, " x ", colBlockSize))
    ## calculate start - and stop position for every block
    if (nrow(data) < rowBlockSize | rowBlockSize == 0) {
        rowBlockSize <- nrow(data)
    }
    if (ncol(data) < colBlockSize | colBlockSize == 0) {
        colBlockSize <- ncol(data)
    }
    rowPos <- calcPositions(nrow(data), rowBlockSize)
    colPos <- calcPositions(ncol(data), colBlockSize)
    blockFrame <- calcBlockFrame(rowPos, colPos, rowBlockSize, colBlockSize)
    
    
    blocks<-apply(blockFrame,1, FUN = function(x, data, total){
        return(list(blockNr=x[1],block=data[x[2]:x[3], x[4]:x[5]], total=total))}, 
        data = data, total=length(blockFrame$number))
    
    
    ## only block numbers needed to distribute the blocks onto
    ## the worker slaves
    blockNumbers <- blockFrame[, 1]
    ## run BEclear in parallel mode
    blocksDone <-
        unlist(bplapply(blocks, imputeMissingDataForBlock, dir = dir, epochs = epochs,
                        BPPARAM = BPPARAM, lambda = lambda, gamma = gamma,
                        r = r, matrixOfOnes = matrixOfOnes))
    
    ## combine the blocks to the predictedGenes data.frame
    predictedGenes <- combineBlocks(blockFrame, rowPos, colPos, dir)
    colnames(predictedGenes) <- colnames(data)
    
    ## remove all stored single block files
    blockFilenames <- c()
    for (i in seq_len(nrow(blockFrame))) {
        row <- paste(
            "D",
            blockFrame$number[i],
            ".RData",
            sep = ""
        )
        blockFilenames <- c(blockFilenames, row)
    }
    for (i in seq_len(length(blockFilenames))) {
        filedir <- paste(dir, blockFilenames[i], sep = "/")
        file.remove(filedir)
    }
    
    remove(blockFrame, blockNumbers, colPos, rowPos)
        
    ## save block as predictedGenes
    predictedGenes <- as.data.frame(predictedGenes)
    if (outputFormat == "RData") {
        
        filename = paste("predicted.genes", "RData", sep = ".")
        save(predictedGenes, file = paste(dir, filename, sep = "/"))
        
    } else if (outputFormat == "txt") {
        
        filename = paste("predicted.genes", "txt", sep = ".")
        write.table(predictedGenes, file = paste(dir, filename, sep = "/"),
                    row.names = TRUE, col.names = TRUE, sep = "\t")
        
    }
    return(as.matrix(predictedGenes))
    
}
