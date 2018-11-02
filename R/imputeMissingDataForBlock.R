#' imputeMissingDataForBlock
#' 
#' @import futile.logger
#' @import Matrix
#' @importFrom stats rnorm
#' 
#' @keywords internal
#' 
#' @return number of the block processed
imputeMissingDataForBlock <- function(data, block, blockFrame, dir, epochs, 
                                      lambda = 1, gamma = 0.01, r = 10) {
    
    flog.info(paste("Impute missind data for block", block, "of", 
                    length(blockFrame[[2]])))
    
    rowStartPosition <- blockFrame[[2]][[block]]
    rowStopPosition <- blockFrame[[3]][[block]]
    colStartPosition <- blockFrame[[4]][[block]]
    colStopPosition <- blockFrame[[5]][[block]]
    
    
    ## take one block of data
    D<-data[rowStartPosition:rowStopPosition,
            colStartPosition:colStopPosition]
    
    ## check if NA values are contained in the block
    if (any(is.na(D))) {
        
        
        ## set D as matrix
        dat <- as.matrix(D)
        
        ## set NA data to 0
        D <- dat
        D[is.na(D)] <- 0
        D <- as(as.matrix(D), "dgCMatrix")
        
        Dsummary <- summary(D)
        m <- nrow(D)               # number of rows
        n <- ncol(D)               # number of columns
        is <- Dsummary$i           # row of each revealed entry
        js <- Dsummary$j           # column of each revealed entry
        xs <- Dsummary$x           # value of each revealed entry
        N <- length(is)            # number of revealed entries
        
        set.seed(1, kind="Mersenne-Twister")
        L0r10 <- matrix(rnorm(m * r), m, r) / sqrt(r)
        R0r10 <- matrix(rnorm(r * n), r, n) / sqrt(r)
        ## run LFM 
        resultGdr10l1 <- runGradientDescent(L0r10, R0r10, lambda = lambda, 
                                            epochs = epochs, gamma = gamma, 
                                            block = block, is = is, js = js, 
                                            D = dat, r = r)
        
        dataTemp <- data[rowStartPosition:rowStopPosition,
                         colStartPosition:colStopPosition]
        ## Predicted matrix
        
        D1 <- resultGdr10l1$L %*% resultGdr10l1$R
        ## Replace predicted values with known, keep predicted unknown
        ## values
        for(i1 in seq_len(nrow(dataTemp)))
            for(j1 in seq_len(ncol(dataTemp)))
            {
                if(!is.na(dataTemp[i1, j1]))
                {
                    D1[i1, j1] <- dataTemp[i1, j1]
                }
            }
        colnames(D1) <- colnames(dataTemp)
        rownames(D1) <- rownames(dataTemp)
        
        blockName <- paste("D", block, sep="")
        
        filename <- paste(blockName, "row", rowStartPosition,
                          rowStopPosition, "col", colStartPosition, 
                          colStopPosition,"RData", sep=".")
        save(D1, file=paste(dir, filename, sep="/"))
        
        remove(D1, D, dataTemp)
        return(block)
    }
    
    ## no NA values contained in the block - keep original values
    else {
        flog.debug(paste("Block", block, 
                         "has no missing values. Original values are kept"))
        D1 <- D
        colnames(D1) <- colnames(D)
        rownames(D1) <- rownames(D)
        
        blockName <- paste("D", block, sep="")
        filename <- paste(blockName, "row", rowStartPosition, 
                          rowStopPosition, "col", colStartPosition, 
                          colStopPosition, "RData", sep=".")
        save(D1, file=paste(dir, filename, sep="/"))
        
        remove(D1, D)
        return(block)
    }
}