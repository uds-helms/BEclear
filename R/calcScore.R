#' calcScore
#' 
#' @aliases calcScore
#' 
#' @title calculate batch effect score
#' 
#' @seealso \code{\link{calcMedians}}
#' @seealso \code{\link{calcPvalues}}
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
#' @param summary a summary data.frame containing the columns "gene", "batch", 
#' "median" and "p-value" and consists of all genes which were found in the 
#' median and p-value calculations, see \code{\link{calcSummary}} function for 
#' more details.
#' @param dir set the path to a directory the returned data.frame should be 
#' stored. The current working directory is defined as default parameter.
#' 
#' @export calcScore
#' @import futile.logger
#' @usage calcScore(data, samples, summary, dir=getwd())
#' 
#' @return A data.frame is returned containing the number of found genes assumed
#' to be batch affected separated by batch and a BEscore for every batch. The 
#' data.frame is also stored in the specified directory as .RData file.
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
#' # Calculates median difference values and p-values from the example data
#'library(data.table)
#'samples <- data.table(ex.samples)
#'data <- data.table(feature=rownames(ex.data), ex.data)
#'data <- melt(data = data, id.vars = "feature", variable.name = "sample", 
#'             value.name = "beta.value")
#'setkey(data, "feature", "sample")
#'med <- calcMedians(data=data, samples=samples)
#' pvals <- calcPvalues(data=ex.data, samples=ex.samples, adjusted=TRUE, 
#' method="fdr")
#' 
#' # Summarize p-values and median differences for batch affected genes 
#' sum <- calcSummary(medians=med, pvalues=pvals)
#' 
#' # Calculates the score table
#' score.table <- calcScore(data=ex.data, samples=ex.samples, summary=sum,
#' dir=getwd())

calcScore <- function(data, samples, summary, dir=getwd()) {
    ## take batch numbers
    batches <- unique(samples$batch_id)
    flog.info(paste("Calculating the scores for", length(batches), "batches"))
    ## take number of genes
    numGenes <- nrow(data)
    ## table with number of found genes & score for every batch
    geneTableMedians <- as.data.frame(matrix(ncol=12))
    colnames(geneTableMedians) <- c("batch", "0.05", "0.1", "0.2", "0.3",
                                    "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", 
                                    "BEscore")
    counter <- count005 <- count01 <- count02 <- count03 <- count04 <-
        count05 <- count06 <- count07 <- count08 <- count09 <- 0
    for (s in batches) {
        counter <- counter + 1
        for (i in seq_len(nrow(summary))) {
            ## values belong to batch s
            if (summary[i, "batch"] == s & summary[i, "pvalue"] <= 0.01) {
                
                if ((is.na(summary[i, "median"]) 
                     | summary[i, "median"] >= 0.9)) {
                    count09 <- count09 + 1
                }   else if (summary[i, "median"] >= 0.8) {
                    count08 <- count08 + 1
                }   else if (summary[i, "median"] >= 0.7) {
                    count07 <- count07 + 1
                }   else if (summary[i, "median"] >= 0.6) {
                    count06 <- count06 + 1
                }   else if (summary[i, "median"] >= 0.5) {
                    count05 <- count05 + 1
                }   else if (summary[i, "median"] >= 0.4 ) {
                    count04 <- count04 + 1
                }   else if (summary[i, "median"] >= 0.3) {
                    count03 <- count03 + 1
                }   else if (summary[i, "median"] >= 0.2) {
                    count02 <- count02 + 1
                }   else if (summary[i, "median"] >= 0.1) {
                    count01 <- count01 + 1
                }   else if (summary[i, "median"] >= 0.05) {
                    count005 <- count005 + 1
                } 
            }
        }
        geneTableMedians[counter, "batch"] <- s
        geneTableMedians[counter, "0.05"] <- count005
        geneTableMedians[counter, "0.1"] <- count01
        geneTableMedians[counter, "0.2"] <- count02
        geneTableMedians[counter, "0.3"] <- count03
        geneTableMedians[counter, "0.4"] <- count04
        geneTableMedians[counter, "0.5"] <- count05
        geneTableMedians[counter, "0.6"] <- count06
        geneTableMedians[counter, "0.7"] <- count07
        geneTableMedians[counter, "0.8"] <- count08
        geneTableMedians[counter, "0.9"] <- count09
        beScore <- count005*1 + count01*2 + count02*4 + count03*6 + 
            count04*8 + count05*10 + count06*12 + count07*14 + count08*16 + 
            count09*18
        if (numGenes > 0) {
            beScore <- beScore / numGenes
        }
        geneTableMedians[counter, "BEscore"] <- beScore
        count005 <- count01 <- count02 <- count03 <- count04 <- count05 <- 
            count06 <- count07 <- count08 <- count09 <- 0
    }
    
    remove(counter, count005, count01, count02, count03, count04, count05, 
           count06, count07, count08, count09, beScore, i, s)
    scoreTable <- geneTableMedians
    remove(geneTableMedians)
    filename="score.table.Rdata"
    save(scoreTable, file=paste(dir, filename, sep="/"))
    return(scoreTable)
}
