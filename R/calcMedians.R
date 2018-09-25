#' calcMedians
#' 
#' @aliases calcMedians
#' 
#' @title Calculate median difference values
#' 
#' @seealso \code{\link{correctBatchEffect}}
#' 
#' @description Compares the median value of all beta values belonging to one 
#' batch with the median value of all beta values belonging to all other 
#' batches. Returns a matrix containing this median difference value for every 
#' gene in every batch, columns define the batch numbers, rows the gene names.
#' 
#' 
#' @param data any matrix filled with beta values, column names have to be 
#' sample_ids corresponding to the ids listed in "samples", row names have to be
#' gene names.
#' @param samples data frame with two columns, the first column has to contain 
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param BPPARAM An instance of the 
#' \code{\link[BiocParallel]{BiocParallelParam}} that determines how to 
#' parallelisation of the functions will be evaluated.
#' 
#' @export calcMedians
#' @import BiocParallel
#' @import futile.logger
#' @usage calcMedians(data, samples, BPPARAM=bpparam())
#' 
#' @return a matrix containing median comparison values for all genes in all 
#' batches, the column names define the batch numbers, row names are the same 
#' gene names as contained in the input matrix.
#' 
#' @examples 
#' ## Shortly running example. For a more realistic example that takes
#' ## some more time, run the same procedure with the full BEclearData
#' ## dataset.
#' 
#' ## Calculate median comparison values in non-parallel mode
#' data(BEclearData)
#' ex.data <- ex.data[31:90,7:26]
#' ex.samples <- ex.samples[7:26,]
#' medians <- calcMedians(data=ex.data, samples=ex.samples)

calcMedians <- function(data, samples, BPPARAM=bpparam()) {
    ## get batch numbers
    batches <- unique(samples$batch_id)
    flog.info(paste("Calculating the medians for", length(batches), "batches"))
    ## get genes
    genes <- rownames(data)
    
    ## construct data.frames filled with NA, one row per gene, one column per 
    ## batch
    medianDif <- as.data.frame(matrix(NA, nrow=length(genes),
                                      ncol=length(batches)))
    rownames(medianDif) <- genes
    colnames(medianDif) <- batches
    
    result <- bplapply(batches, calcMediansForBatch, genes = genes, 
                       medianDif = medianDif, samples = samples,
                       data = data, BPPARAM=BPPARAM)
    
    ## fill median matrix from result
    result <- unlist(result)
    counter <- 1
    for(j in seq_len(ncol(medianDif))) {
        for(i in seq_len(nrow(medianDif))) {
            medianDif[i, j] <- result[counter]
            counter <- counter + 1
        }
    }
    remove(result, counter, i, j)
    
    return(medianDif)
}
