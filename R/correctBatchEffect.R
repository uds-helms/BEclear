#' correctBatchEffect
#'
#' @aliases correctBatchEffect
#'
#' @seealso \code{\link{calcMedians}}
#' @seealso \code{\link{calcPvalues}}
#' @seealso \code{\link{calcSummary}}
#' @seealso \code{\link{calcScore}}
#' @seealso \code{\link{clearBEgenes}}
#' @seealso \code{\link{imputeMissingData}}
#' @seealso \code{\link{replaceWrongValues}}
#'
#' @title Correct a batch effect in DNA methylation data
#'
#' @description This method combines most functions of the
#' \code{\link{BEclear-package}} to one. The method performs the whole process
#' of searching for batch effects and automatically correct them for a matrix
#' of beta values stemming from DNA methylation data.
#'
#' @details The function performs the whole process of searching for batch
#' effects and automatically correct them for a matrix of beta values stemming
#' from DNA methylation data. Thereby, the function is running most of the
#' functions from the \code{\link{BEclear-package}} in a logical order.\cr
#' First, median comparison values are calculated by the
#' \code{\link{calcMedians}} function, followed by the calculation of p-values
#' by the \code{\link{calcPvalues}} function. With the results from the median
#' comparison and p-value calculation, a summary data frame is build using the
#' \code{\link{calcSummary}} function, and a scoring table is established by
#' the \code{\link{calcScore}} function. Now, found entries from the summary are
#' set to NA in the input matrix using the \code{\link{clearBEgenes}} function,
#' then the \code{\link{imputeMissingData}} function is used to predict the
#' missing values and at the end, possibly existing wrongly predicted entries
#' (values lower than 0 or greater than 1) are corrected using the
#' \code{\link{replaceWrongValues}} function.
#'
#' @references \insertRef{Akulenko2016}{BEclear}
#' @references \insertRef{Koren2009}{BEclear}
#' @references \insertRef{Candes2009}{BEclear}
#'
#' @param data any matrix filled with beta values, column names have to be
#' sample_ids corresponding to the ids listed in "samples", row names have to
#' be gene names.
#' @param samples data frame with two columns, the first column has to contain
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param adjusted should the p-values be adjusted or not, see "method" for
#' available adjustment methods.
#' @param method adjustment method for p-value adjustment (if TRUE), default
#' method is "false discovery rate adjustment", for other available methods see
#'  the description of the used standard R package \code{\link{p.adjust}}. See
#' \code{\link{calcPvalues}} for more information
#' @param rowBlockSize the number of rows that is used in a block if the
#' function is run in parallel mode and/or not on the whole matrix. Set this,
#' and the "colBlockSize" parameter to 0 if you want to run the function on the
#' whole input matrix. See \code{\link{imputeMissingData}} and especially the
#' details section of See \code{\link{imputeMissingData}} for more information
#' about this feature.
#' @param colBlockSize the number of columns that is used in a block if the
#' function is run in parallel mode and/or not on the whole matrix. Set this,
#' and the "rowBlockSize" parameter to 0 if you want to run the function on the
#' whole input matrix. See \code{\link{imputeMissingData}} and especially the
#' details section of See \code{\link{imputeMissingData}} for more information
#' about this feature.
#' @param epochs the number of iterations used in the gradient descent algorithm
#' to predict the missing entries in the data matrix. See 
#' \code{\link{imputeMissingData}} for more information.
#' @param lambda constant that controls the extent of regularization during the 
#' gradient descent
#' @param gamma constant that controls the extent of the shift of parameters 
#' during the gradient descent
#' @param r length of the second dimension of variable matrices R and L
#' @param outputFormat you can choose if the finally returned data matrix should
#' be saved as an .RData file or as a tab-delimited .txt file in the specified
#' directory. Allowed values are "RData" and "txt".
#' See \code{\link{imputeMissingData}}
#' for more information.
#' @param dir set the path to a directory the predicted matrix should be stored.
#' The current working directory is defined as default parameter.
#' @param BPPARAM An instance of the
#' \code{\link[BiocParallel]{BiocParallelParam-class}} that determines how to
#' parallelisation of the functions will be evaluated.
#'
#' @export correctBatchEffect
#' @import BiocParallel
#' @import futile.logger
#' @import data.table
#' @usage correctBatchEffect(data, samples, adjusted=TRUE, method="fdr",
#' rowBlockSize=60, colBlockSize=60, epochs=50, lambda = 1, gamma = 0.01, r = 10,
#' outputFormat="RData", dir=getwd(), BPPARAM=bpparam())
#'
#' @return A list containing the following fields (for detailed information look
#' at the documentations of the corresponding functions):
#' \describe{
#' \item{medians}{A data.frame containing all median comparison values
#' corresponding to the input matrix.}
#' \item{pvalues}{A data.frame containing all p-values corresponding to the
#' input matrix.}
#' \item{summary}{The summarized results of the median comparison and p-value
#' calculation.}
#' \item{score.table}{A data.frame containing the number of found genes and a
#' BEscore for every batch.}
#' \item{cleared.data}{the input matrix with all values defined in the summary
#' set to NA.}
#' \item{predicted.data}{the input matrix after all previously NA values have
#' been predicted.}
#' \item{corrected.predicted.data}{the predicted matrix after the correction for
#' wrongly predicted values.}
#' }
#'
#' @examples
#' ## Shortly running example. For a more realistic example that takes
#' ## some more time, run the same procedure with the full BEclearData
#' ## dataset.
#'
#' ## Whole procedure that has to be done to use this function.
#' ## Correct the example data for a batch effect
#' data(BEclearData)
#' ex.data <- ex.data[31:90,7:26]
#' ex.samples <- ex.samples[7:26,]
#'
#' # Note that row- and block sizes are just set to 10 to get a short runtime.
#' # To use these parameters, either use the default values or please note the
#' # description in the details section of \code{\link{imputeMissingData}}
#' result <- correctBatchEffect(data=ex.data, samples=ex.samples,
#' adjusted=TRUE, method="fdr", rowBlockSize=10, colBlockSize=10,
#' epochs=50, outputFormat="RData", dir=getwd())
#'
#' # Unlist variables
#' medians <- result$medians
#' pvals <- result$pvals
#' summary <- result$summary
#' score <- result$score.table
#' cleared <- result$clearedData
#' predicted <- result$predictedData
#' corrected <- result$correctedPredictedData

correctBatchEffect <- function(data, samples, adjusted = TRUE, method = "fdr",
                               rowBlockSize = 60, colBlockSize = 60, 
                               epochs = 50, lambda = 1, gamma = 0.01, r = 10,
                               outputFormat = "RData",
                               dir = getwd(), BPPARAM = bpparam()) {
    
    ## checking if they're are values above 1
    if(any(data > 1, na.rm=TRUE)){
        flog.warn(paste(sum(data > 1, na.rm = TRUE), "values are above 1"))
        flog.warn("Replacing them with missing values")
        data[data > 1] <- NA
    }
    
    ## checking if there are columns containing only missing values
    naIndices <- apply(data, 2, function(x) all(is.na(x)))
    if(any(naIndices, na.rm=TRUE)){
        flog.warn("There are columns, that contain only missing values")
        flog.warn(paste(sum(naIndices), "columns get dropped"))
        data <- data[ ,!naIndices]
    }
    
    ## checking if there are rows containing only missing values
    naIndices <- apply(data, 1, function(x) all(is.na(x)))
    if(any(naIndices, na.rm=TRUE)){
        flog.warn("There are rows, that contain only missing values")
        flog.warn(paste(sum(naIndices), "rows get dropped"))
        data <- data[!naIndices, ]
    }
    
    ## checking if there are samples that are not present in the samples matrix
    if(any(!colnames(data) %in% samples$sample_id)){
        ids <- paste(colnames(data)[!colnames(data) %in% samples$sample_id], 
                     collapse= ", ")
        flog.warn(paste("The following samples are in the data, but not annotated", 
                        "in the samples matrix:", ids))
        flog.warn("Dropping those samples now")
        data <- data[,colnames(data) %in% samples$sample_id]
    }
    
    if(any(!samples$sample_id %in% colnames(data))){
        ids <- paste(samples$sample_id[!samples$sample_id %in% colnames(data)],
                     collapse = ", ")
        flog.warn("The following samples are annotated in the sample matrix,",
                  "but aren't contained in data matrix:", ids)
        flog.warn("Dropping those samples now")
        samples <- samples[sample_id %in% colnames(data)]
    }
    
    samples <- data.table(samples)
    uniqueIDsToSamples <- NULL
    
    ##checking if there are duplicated sample names
    if(any(duplicated(colnames(data)))){
        flog.warn("Sample names aren't unique")
        flog.warn(paste("Transforming them to unique IDs. List with annotations will",
                  "be added to the results"))
        uniqueIDsToSamples <- data.table(sample_id = colnames(data), 
                                         unique_id = 
                                             as.character(seq_along(colnames(data))))
        colnames(data) <- uniqueIDsToSamples$unique_id
        samples <- samples[uniqueIDsToSamples, ,
                           on=.(sample_id = sample_id)][,
                                                        .(sample_id = unique_id,
                                              batch_id = batch_id)]
    }
    setkey(samples, "batch_id", "sample_id")
    
    flog.info("Transforming matrix to data.table")
    DT <- data.table(feature=as.character(rownames(data)), data)
    DT <- melt(data = DT, id.vars = "feature", variable.name = "sample", 
                 value.name = "beta.value", variable.factor = FALSE)
    setkey(DT, "feature", "sample")
    
    med <- calcMedians(DT, samples, BPPARAM = BPPARAM)
    pval <- calcPvalues(DT, samples, adjusted, method, BPPARAM = BPPARAM)
    sum <- calcSummary(med, pval)
    
    #flog.info("Transforming data.table back to matrix")
    #data <- dcast(data, feature ~ sample, value.var = "beta.value")
    #data <- as.matrix(data, rownames = "feature")
    rm(DT)
    if(is.null(sum)){
        flog.info("There were no batch effects detected")
        score <- NULL
        cleared <- data
    }else{
        score <- calcScore(data, samples, sum)
        cleared <- clearBEgenes(data, samples, sum)
    }
    
    predicted <-
        imputeMissingData (data=cleared, rowBlockSize=rowBlockSize, 
                           colBlockSize=colBlockSize, epochs=epochs,
                           lambda =lambda, gamma = gamma, r = r,
                           outputFormat = outputFormat, dir=dir, BPPARAM = BPPARAM)
    corrected <- replaceWrongValues(predicted)
    
    return(list(medians = med, pvals = pval, summary = sum, 
                scoreTable = score, clearedData = cleared, 
                predictedData = predicted, correctedPredictedData = 
                    corrected, uniqueIDsToSamples = uniqueIDsToSamples))
}