#' calcBatchEffects
#'
#' @aliases calcBatchEffects
#'
#' @title Calculate the Batch Effects in a given data set
#'
#' @seealso \code{\link{ks.test}}
#' @seealso \code{\link{p.adjust}}
#' @seealso \code{\link{correctBatchEffect}}
#'
#' @description Calculates for each gene in every batch the median distance to the other
#' batches and the p-value resulting from the Kolmogorov-Smirnov test.
#'
#' @details
#' \enumerate{
#' \item{medians}
#' Compares the median value of all beta values belonging to one
#' batch with the median value of all beta values belonging to all other
#' batches. Returns a matrix containing this median difference value for every
#' gene in every batch, columns define the batch numbers, rows the gene names.
#' \item{p-values}
#' Compares the distribution of all beta values corresponding to
#' one batch with the distribution of all beta values corresponding to all
#' other batches and returns a p-value which defines if the distributions are
#' the same or not. Standard two sided Kolmogorov-Smirnov test is used to
#' calculate the (adjusted) p-values.
#' }
#'
#' @param data a \code{\link[data.table]{data.table}} with one column indicating
#' the sample, one the features and a value column indicating the beta value
#' or a matrix with rows as features and columns as samples
#' @param samples data frame with two columns, the first column has to contain
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param adjusted should the p-values be adjusted or not, see "method" for
#' available adjustment methods.
#' @param method adjustment method for p-value adjustment (if TRUE), default
#' method is "false discovery rate adjustment", for other available methods
#' see the description of the used standard R package \code{\link{p.adjust}}.
#' @param BPPARAM An instance of the
#' \code{\link[BiocParallel]{BiocParallelParam-class}} that determines how to
#' parallelisation of the functions will be evaluated.
#'
#' @export calcBatchEffects
#' @import futile.logger
#' @import data.table
#' @importFrom stats p.adjust
#' @importFrom methods is
#' @usage calcBatchEffects(data, samples, adjusted=TRUE, method="fdr",
#' BPPARAM=SerialParam())
#'
#' @return a matrix containing medians and p-values for all genes in all batches
#' @examples
#' ## Shortly running example. For a more realistic example that takes
#' ## some more time, run the same procedure with the full BEclearData
#' ## dataset.
#' 
#' ## Calculate fdr-adjusted p-values in non-parallel mode
#' data(BEclearData)
#' ex.data <- ex.data[31:90, 7:26]
#' ex.samples <- ex.samples[7:26, ]
#' 
#' 
#' res <- calcBatchEffects(data = ex.data, samples = ex.samples, method = "fdr")
calcBatchEffects <- function(data, samples, adjusted = TRUE, method = "fdr",
                             BPPARAM = SerialParam()) {
  if (!is(data, "data.table")) {
    flog.info("Transforming matrix to data.table")
    data <- data.table(feature = as.character(rownames(data)), data)
    data <- melt(
      data = data, id.vars = "feature", variable.name = "sample",
      value.name = "beta.value", variable.factor = FALSE
    )
    setkey(data, "feature", "sample")
  }

  if (!is(samples, "data.table")) {
    samples <- data.table(samples)
  }

  flog.info(paste(
    "Calculating the batch effects for",
    samples[, uniqueN(batch_id)], "batches"
  ))

  batchEffects <- lapply(samples[, unique(batch_id)], calcBatchEffectsForBatch,
    samples = samples, data = data, BPPARAM = BPPARAM)

  flog.debug(paste("Binding", length(batchEffects), "columns of batch_effects together"))
  batchEffects <- do.call(cbind, batchEffects)

  med <- batchEffects[colnames(batchEffects) == "medians"]
  colnames(med) <- samples[, unique(batch_id)]

  pvalues <- batchEffects[colnames(batchEffects) == "p_values"]
  colnames(pvalues) <- samples[, unique(batch_id)]

  # pvalue adjustment
  if (adjusted == TRUE) {
    flog.info("Adjusting p-values")
    pvalues <- t(apply(pvalues, 1, p.adjust, method = method))
  }

  return(list(med = med, pval = pvalues))
}
