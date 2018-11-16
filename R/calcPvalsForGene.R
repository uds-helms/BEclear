#' calcPvalsForGene
#'
#'@keywords internal
#'
#'@description  function to calculate p-values for one gene in a batch
#'
#'@importFrom stats ks.test
#'@import futile.logger
#'
#'@return the p-values for genes in a batch
calcPvalsForGene <- function(data){
    
    if(all(is.na(data$batch))){
        return(0.0)
    }else{
        return(ks.test(data$batch, data$others)$p.value)
    }
}