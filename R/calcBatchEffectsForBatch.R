#' calcBatchEffectsForBatch
#'
#' @keywords internal
#'
#' @description  function to calculate batch effects for every gene in a batch
#'
#' @importFrom stats ks.test
#' @import BiocParallel
#' @import futile.logger
#'
#' @return the medians p-values for genes in a batch
calcBatchEffectsForBatch <- function(batch, samples, data, 
                                     BPPARAM = SerialParam()) {
  flog.debug(paste("Calculating the batch effect for batch", batch))

  flog.debug("Subsetting data.table for batch")
  DT_batch <- samples[batch_id == batch][data
    , ,
    on = .(sample_id = sample),
    nomatch = 0
  ][, .(feature, beta.value)]
  flog.debug("Subsetting data.table for others")
  DT_other <- samples[batch_id != batch][data
    , ,
    on = .(sample_id = sample),
    nomatch = 0
  ][, .(feature, beta.value)]


  features_batch <- DT_batch[, list(list(beta.value)), by = feature]$V1
  features_other <- DT_other[, list(list(beta.value)), by = feature]$V1

  flog.debug(paste("Calculating the p-values for batch", batch))
  result <- bpmapply(function(X, Y) {
    medianDiff <- median(X, na.rm = TRUE) - median(Y, na.rm = TRUE)
    if (all(is.na(X))) {
      return(c(0.0, medianDiff))
    } else {
      suppressWarnings(pval <- ks.test(X, Y)$p.value)
      return(c(pval, medianDiff))
    }
  }, X = features_batch, Y = features_other, BPPARAM = BPPARAM)

  result <- t(result)
  colnames(result) <- c("p_values", "medians")

  return(data.frame(result, row.names = unique(data$feature)))
}
