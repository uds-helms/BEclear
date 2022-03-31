#' calcScore
#'
#' @aliases calcScore
#'
#' @title calculate batch effect score
#'
#' @seealso \code{\link{calcBatchEffects}}
#' @seealso \code{\link{calcSummary}}
#' @seealso \code{\link{correctBatchEffect}}
#'
#' @description Returns a table with the number of found genes with found
#' p-values less or equal to 0.01 and median values greater or equal to 0.05.
#' A score is calculated depending on the number of found genes as well as the
#' magnitude of the median difference values, this score is divided by the
#' overall number of genes in the data and returned as "BEscore". See details
#' for further information and details about the score calculation. The
#' returned data.frame is also stored in the specified directory as .RData file.
#'
#' @details The returned data frame contains one column for the batch numbers,
#' 11 columns containing the number of genes found in a certain range of the
#' median difference value and a column with the calculated BEscore. These
#' found genes are assumed to be batch affected due to their difference in
#' median values and their different distribution of the beta values. The higher
#' the found number of genes and the more extreme the median difference is, the
#' more severe is the assumed batch effect supposed to be. We suggest that there
#' is no need for a batch effect correction if the BEscore for a batch is less
#' than 0.02. BEscores between 0.02 and 0.1 are lying in a "gray area" for which
#' we assume a not severe batch effect, and values beyond 0.1 certainly describe
#' a batch effect and should definitely be corrected.\cr
#' The 11 columns containing the numbers of found genes count the median
#' difference values which are ranging from >= 0.05 to < 0.1 ; >= 0.1 to < 0.2;
#' >= 0.2 to < 0.3 and so on up to a limit of 1.\cr
#' The BEscore is calculated by the sum of the weighted number of genes divided
#' by the number of genes. Weightings are calculated by multiplication of the
#' number of found genes between 0.05 and 0.1 by 1, between 0.1 and 0.2 by 2,
#' between 0.2 and 0.3 by 4, between 0.3 and 0.4 by 6 and so on.
#'
#' @param data any matrix filled with beta values, column names have to be
#' sample_ids corresponding to the ids listed in "samples", row names have to
#' be gene names.
#' @param samples data frame with two columns, the first column has to contain
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param summary a summary \code{\link[data.table]{data.table}} containing the 
#' columns "gene", "batch", median" and "p-value" and consists of all genes which 
#' were found in the median and p-value calculations, see 
#' \code{\link{calcSummary}} function for more details.
#' @param saveAsFile determining if the data.frame should also be saved as a file
#' @param dir set the path to a directory the returned data.frame should be
#' stored. The current working directory is defined as default parameter.
#'
#' @export calcScore
#' @import futile.logger
#' @import data.table
#' @importFrom dixonTest dixonTest
#' @usage calcScore(data, samples, summary, saveAsFile=FALSE, dir=getwd())
#' 
#' @references \insertRef{Dixon1950}{BEclear}
#' @references \insertRef{Dixon1951}{BEclear}
#' @references \insertRef{Rorabacher1991}{BEclear}
#'
#' @return A data.frame is returned containing the number of found genes assumed
#' to be batch affected separated by batch and a BEscore for every batch. Furthermore
#' there's a column dixonPval giving you a p-value regarding each BEscore according 
#' to a Dixon test.
#' The data.frame is also stored in the specified directory as .RData file, if 
#' saveAsFile is \code{TRUE}.
#'
#' @examples
#' ## Shortly running example. For a more realistic example that takes
#' ## some more time, run the same procedure with the full BEclearData
#' ## dataset.
#' 
#' ## Whole procedure that has to be done to use this function.
#' data(BEclearData)
#' ex.data <- ex.data[31:90, 7:26]
#' ex.samples <- ex.samples[7:26, ]
#' ## Calculate the batch effects
#' batchEffects <- calcBatchEffects(data = ex.data, samples = ex.samples,
#' adjusted = TRUE, method = "fdr")
#' med <- batchEffects$med
#' pvals <- batchEffects$pval
#' 
#' # Summarize p-values and median differences for batch affected genes
#' sum <- calcSummary(medians = med, pvalues = pvals)
#' 
#' # Calculates the score table
#' score.table <- calcScore(data = ex.data, samples = ex.samples, summary = sum)
calcScore <- function(data, samples, summary, saveAsFile = FALSE, dir = getwd()) {
  ## take batch ids
  batches <- unique(samples$batch_id)
  flog.info(paste("Calculating the scores for", length(batches), "batches"))
  
  ## take number of genes
  numGenes <- nrow(data)
  
  ## count the medians
  DT <- lapply(batches, FUN = function(x, y){
      data.table(batch_id = x, 
                 count05 = y[batch_id == x & median >= 0.05 & median < 0.1, .N], 
                 count1 = y[batch_id == x & median  >= 0.1  & median < 0.2,  .N], 
                 count2 = y[batch_id == x & median  >= 0.2  & median < 0.3,  .N], 
                 count3 = y[batch_id == x & median  >= 0.3  & median < 0.4,  .N],
                 count4 = y[batch_id == x & median  >= 0.4  & median < 0.5,  .N], 
                 count5 = y[batch_id == x & median  >= 0.5  & median < 0.6,  .N], 
                 count6 = y[batch_id == x & median  >= 0.6  & median < 0.7,  .N], 
                 count7 = y[batch_id == x & median  >= 0.7  & median < 0.8,  .N], 
                 count8 = y[batch_id == x & median  >= 0.8  & median < 0.9,  .N], 
                 count9 = y[batch_id == x & median  >= 0.9,  .N])}, summary)
  
  DT <- do.call(rbind, DT)
  
  ## calculate BEscores
  DT[, BEscore := (count05 * 1 + count1* 2 + count2 * 4 + count3 * 6 + count4 * 8 
                   + count5 * 10 + count6 * 12 + count7 * 14 + count8 * 16 
                   + count9 * 18)/numGenes ]
  
  
  
  ## calculate outlier according to the dixon test
  DT[, dixonPval := as.numeric(NA)]
  if(sum(!is.na(DT$BEscore)) >= 3 && max(DT$BEscore) != min(DT$BEscore)){
    pval <- dixonTest(DT$BEscore, alternative = c("two.sided"))
    DT[BEscore == max(DT$BEscore), dixonPval := as.numeric(pval$p.value)]
  }
  
  
  
  if (saveAsFile) {
      filename <- "score.table.Rdata"
      save(DT, file = paste(dir, filename, sep = "/"))
  }
  
  return(DT)
}
