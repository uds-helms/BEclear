#' BEclear-package
#'
#' @useDynLib BEclear
#'
#' @aliases BEclear-package
#' @aliases BEclearCorrected
#'
#' @docType package
#'
#' @title Correction of batch effects in DNA methylation data
#'
#'
#' @description Provides some functions to detect and correct for batch effects
#' in DNA methylation data. The core function \code{\link{correctBatchEffect}} is based on
#' Latent Factor Models and can also be used to predict missing values in any
#' other matrix containing real numbers.
#'
#' @details \code{\link{correctBatchEffect}}:
#' The function combines most functions of the \code{\link{BEclear-package}} to
#' one. This function performs the whole process of searching for batch effects
#' and automatically correct them for a matrix of beta values stemming from DNA
#' methylation data.\cr
#' \code{\link{correctBatchEffect}}:
#' This function predicts the missing entries of an input matrix (NA values)
#' through the use of a Latent Factor Model.\cr
#' \code{\link{calcBatchEffects}}:
#' Compares the median value of all beta values belonging to one batch with the
#' median value of all beta values belonging to all other batches. Returns a
#' matrix containing this median difference value for every gene in every batch,
#' columns define the batch numbers, rows the gene names.\cr
#' And compares the distribution of all beta values corresponding to one batch 
#' with the distribution of all beta values corresponding to all other batches and
#' returns a p-value which defines if the distributions are the same or not.\cr
#' \code{\link{calcSummary}}:
#' Summarizes the results of the \code{\link{calcBatchEffects}} function\cr
#' \code{\link{calcScore}}:
#' Returns a table with the number of found genes with found p-values less or
#' equal to 0.01 and median values greater or equal to 0.05. A score is
#' calculated depending on the number of found genes as well as the magnitude of
#' the median difference values, this score is divided by the overall number of
#' genes in the data and returned as "BEscore". See the methods details for
#' further information and details about the score calculation.\cr
#' \code{\link{makeBoxplot}}:
#' A simple \code{\link{boxplot}} is done with boxes either separated by batches
#' or by samples and describe the five number summary of all beta values
#' corresponding to a batch or a sample, respectively. The batch_ids are shown on
#' the x-axis with a coloring corresponding to the BEscore.\cr
#' \code{\link{clearBEgenes}}:
#' A function that simply sets all values to NA which were previously found by
#' median value comparison and p-value calculation and are stored in a summary.
#' The summary defines which values in the data matrix are set to NA.\cr
#' \code{\link{countValuesToPredict}}:
#' Simple function that counts all values in a matrix which are NA\cr
#' \code{\link{findOutsideValues}}:
#' A method which lists values below 0 or beyond 1 contained in the input matrix.
#' These entries are stored in a data.frame together with the corresponding
#' row and column position of the matrix.\cr
#' \code{\link{replaceOutsideValues}}:
#' A method which replaces values below 0 or beyond 1 contained in the input
#' matrix. These entries outside the boundaries are replaced by 0 or 1, respectively.
#'
#' @examples
#' data(BEclearData)
#' ## Calculate the batch effects
#' batchEffects <- calcBatchEffects(data = ex.data, samples = ex.samples,
#' adjusted = TRUE, method = "fdr")
#' med <- batchEffects$med
#' pvals <- batchEffects$pval
#' 
#' ## Summarize p-values and median differences for batch affected genes
#' sum <- calcSummary(medians = med, pvalues = pvals)
#' 
#' ## Calculates the score table
#' score.table <- calcScore(data = ex.data, samples = ex.samples, summary = sum)
#' 
#' ## Simple boxplot for the example data separated by batch
#' makeBoxplot(
#'   data = ex.data, samples = ex.samples, score = score.table,
#'   bySamples = FALSE, main = "Some box plot"
#' )
#' 
#' ## Simple boxplot for the example data separated by samples
#' makeBoxplot(
#'   data = ex.data, samples = ex.samples, score = score.table,
#'   bySamples = TRUE, main = "Some box plot"
#' )
#' 
#' ## Sets assumed batch affected entries to NA
#' cleared <- clearBEgenes(data = ex.data, samples = ex.samples, summary = sum)
#' ## Counts and stores number of entries to predict
#' numberOfEntries <- countValuesToPredict(data = cleared)
#' \dontrun{
#' ## Predicts the missing entries
#' predicted <- imputeMissingData(data = cleared)
#' 
#' ## Find predicted entries outside the boundaries
#' outsideEntries <- findOutsideValues(data = predicted)
#' 
#' ## Replace predicted entries outside the boundaries
#' corrected <- replaceOutsideValues(data = predicted)
#' }
#' 
#' @author Ruslan Akulenko, Markus Merl, Livia Rasp
#'
#' @references \insertRef{Akulenko2016}{BEclear}
#' @import Rdpack
#' @importFrom Rdpack reprompt
"_PACKAGE"
utils::globalVariables(c(
  "batch_id", "beta.value", "feature", "sample_id", ".",
  "unique_id", "dixonPval", "BEscore", "count05", "count1", 
  "count2", "count3", "count4", "count5", "count6", "count7", 
  "count8", "count9"),
package = "BEclear", add = FALSE
)

#' @name BEclear example methylation data
#'
#' @aliases ex.data
#'
#' @docType data
#'
#' @title Example data set for the BEclear-package
#'
#' @usage data(BEclearData)
#'
#' @description Example data set for the BEclear-package
#'
#' @format An example data matrix that is filled with beta values originally
#' stemming from breast cancer data from the TCGA portal [1], colnames are
#' sample ids, rownames are gene names. Generally, beta values are calculated by
#' dividing the methylated signal by the sum of the unmethylated and methylated
#' signals from a DNA methylation microrarray. The sample data used here
#' contains averaged beta values of probes that belong to promoter regions of
#' single genes. Another possibility would be to use beta values of single
#' probes, whereby the probe names should then be used instead of the gene names
#' as rownames of the matrix.
#'
#' @references \insertRef{TCGA}{BEclear}
#'
"ex.data"

#' @name BEclear example sample data
#'
#' @aliases ex.samples
#'
#' @docType data
#'
#' @title Example data set for the BEclear-package
#'
#' @usage data(BEclearData)
#'
#'
#'
#' @format An example data frame containing a column for the sample id and a
#' column for the corresponding batch id, stemming from breast cancer data
#' from the TCGA portal [1]
#'
#' @references \insertRef{TCGA}{BEclear}
"ex.samples"

#' @name ex.corrected.data
#'
#'
#' @docType data
#'
#' @title Example matrix of corrected data for the BEclear-package
#'
#' @description Example matrix containing a already batch effect corrected sample
#' matrix of beta values from breast invasive carcinoma TCGA methylation
#' data.[1] The matrix contains a small amount of predicted beta values outside of
#' the boundaries to show the operating principles of some of the methods from 
#' the BEclear package.
#'
#' @usage data(BEclearCorrected)
#'
#' @format A matrix containing already corrected beta values of some samples from
#' the breast invasive carcinoma TCGA methylation data. The colnames denote
#' samples, rownames denote gene names.
#'
#' @references \insertRef{TCGA}{BEclear}
"ex.corrected.data"
