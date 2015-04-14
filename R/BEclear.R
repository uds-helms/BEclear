calcMedians <- function(data, samples, parallel=TRUE, cores=4) {
    ## get batch numbers
    batches <- unique(samples$batch_id)
    ## get genes
    genes <- rownames(data)

    ## construct data.frames filled with NA, one row per gene, one column per 
    ## batch
    medianDif <- as.data.frame(matrix(NA, nrow=length(genes),
        ncol=length(batches)))
    rownames(medianDif) <- genes
    colnames(medianDif) <- batches

    ## function to calculate medians for every gene in every batch
    calcMedians <- function(batch) {
        dif <- as.data.frame(matrix(NA, nrow=length(genes), ncol=1))
        rownames(dif) <- genes
        for(i in 1:nrow(medianDif)) {
            j <- which(colnames(medianDif) == batch)
            ## barcode ids of all samples from one batch
            batchSamples <- samples$sample_id[samples$batch_id == 
                colnames(medianDif)[j]]
            ## barcode ids of all samples from all other batches
            otherSamples <- samples$sample_id[samples$batch_id != 
                colnames(medianDif)[j]]
            ## save medians for one gene in dif list
            dif[i, 1] <- abs(median(as.numeric(data[rownames(medianDif)[i],
                as.character(batchSamples)])) - 
                median(as.numeric(data[rownames(medianDif)[i],
                as.character(otherSamples)])))
        }
        return(dif)
    }

    if (parallel == TRUE) {
        ## start parallel mode on 4 cpus
        sfInit(parallel=TRUE, cpus=cores, type="SOCK")
        ## load variables needed in parallel computing
        sfExport("batches", "genes", "samples", "data", "medianDif")

        ## calculate medians in parallel mode
        result <- sfClusterApplyLB(batches, calcMedians)

        ## fill median matrix from result
        result <- unlist(result)
        counter <- 1
        for(j in 1:ncol(medianDif)) {
            for(i in 1:nrow(medianDif)) {
                medianDif[i, j] <- result[counter]
                counter <- counter + 1
            }
        }
        remove(result, counter, i, j)

        ## stop parallel mode
        sfStop()
    }
    else {
        ## calc medians for every batch
        for (i in 1:length(batches)) {
            currentBatch <- batches[i]
            result <- calcMedians(currentBatch)
            ## fill median matrix from result
            for(j in 1:nrow(medianDif)) {
                medianDif[j, i] <- result[j, ]
            }
        remove(result, j, currentBatch)
        }
    remove(i)
    }

    return(medianDif)
}

calcPvalues <- function(data, samples, parallel=TRUE, cores=4, adjusted=TRUE,
    method="fdr") {
    ## get batch numbers
    batches <- unique(samples$batch_id)

    ## get genes
    genes <- rownames(data)

    ## construct data.frames filled with NA, one row per gene, one column per 
    ## batch
    pvalues <- as.data.frame(matrix(NA, nrow=length(genes),
        ncol=length(batches)))
    rownames(pvalues) <- genes
    colnames(pvalues) <- batches

    ## function to calculate p-values for every gene in every batch
    calcPvals <- function(batch) {
        pvals <- as.data.frame(matrix(NA, nrow=length(genes), ncol=1))
        rownames(pvals) <- genes
        for(i in 1:nrow(pvalues)) {
            j <- which(colnames(pvalues) == batch)
            ## barcode ids of all samples from one batch
            batchSamples <- samples$sample_id[samples$batch_id == 
                colnames(pvalues)[j]]
            ## barcode ids of all samples from all other batches
            otherSamples <- samples$sample_id[samples$batch_id != 
                colnames(pvalues)[j]]
            ## save pvalue for one gene in pvalue matrix
            pvals[i, 1] <- ks.test(as.numeric(data[rownames(pvalues)[i],
                as.character(batchSamples)]),
                as.numeric(data[rownames(pvalues)[i],
                as.character(otherSamples)]))$p.value
        }
        return(pvals)
    }

    if (parallel == TRUE) {
        ## start parallel mode on 4 cpus
        sfInit(parallel=TRUE, cpus=cores, type="SOCK")
        ## load variables needed in parallel computing
        sfExport("batches", "genes", "samples", "data", "pvalues")

        ## calculate pvalues in parallel mode
        result <- sfClusterApplyLB(batches, calcPvals)

        ## fill pvalue matrix from result
        result <- unlist(result)
        counter <- 1
        for(j in 1:ncol(pvalues)) {
            for(i in 1:nrow(pvalues)) {
                pvalues[i, j] <- result[counter]
                counter <- counter + 1
            }
        }
        remove(result, counter, i, j)

        ## stop parallel mode
        sfStop()
    }
    else {
        ## calc pvalues for every batch
        for (i in 1:length(batches)) {
            currentBatch <- batches[i]
            result <- calcPvals(currentBatch)
            ## fill pvalue matrix from result
            for(j in 1:nrow(pvalues)) {
                pvalues[j, i] <- result[j, ]
            }
            remove(result, j, currentBatch)
        }
    remove(i)
    }

    ## pvalue adjustment
    if (adjusted == TRUE) {
        pvaluesAdjusted <- pvalues
        ## p-value adjustment with false discovery rate
        for(i in 1:nrow(pvaluesAdjusted)) {
            pvaluesAdjusted[i, ] <- p.adjust(pvaluesAdjusted[i, ],
                method=method)
        }
        pvalues <- pvaluesAdjusted
    }

    return(pvalues)
}

calcSummary <- function(medians, pvalues) {
    ## build summary table of found genes
    summaryTable <- as.data.frame(matrix(ncol=4))
    colnames(summaryTable) <- c("gene", "batch", "median", "pvalue")
    counter <-1
    for(i in 1:nrow(medians)) {
        for(j in 1:ncol(medians)) {
            if((medians[i, j] >= 0.05) & (pvalues[i, j] <= 0.01)) {
                summaryTable[counter, "gene"] <- rownames(medians)[i]
                summaryTable[counter, "batch"] <- colnames(medians)[j]
                summaryTable[counter, "median"] <- medians[i, j]
                summaryTable[counter, "pvalue"] <- pvalues[i, j]
                counter <- counter + 1
            }
        }
    }
    return(summaryTable)
}

calcScore <- function(data, samples, summary, dir=getwd()) {
    ## take batch numbers
    batches <- unique(samples$batch_id)
    ## take number of genes
    numGenes <- nrow(data)
    ## table with number of found genes & score for every batch
    geneTableMedians <- as.data.frame(matrix(ncol=12))
    colnames(geneTableMedians) <- c("batch", "0.05", "0.1", "0.2", "0.3",
        "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "BEscore")
    counter <- count005 <- count01 <- count02 <- count03 <- count04 <-
        count05 <- count06 <- count07 <- count08 <- count09 <- 0
    for (s in batches) {
        counter <- counter + 1
        for (i in 1:nrow(summary)) {
            ## values belong to batch s
            if (summary[i, "batch"] == s) {
                if (summary[i, "median"] >= 0.05 & summary[i, "median"] < 0.1 &
                        summary[i, "pvalue"] <= 0.01) {
                    count005 <- count005 + 1
                }
                if (summary[i, "median"] >= 0.1 & summary[i, "median"] < 0.2 &
                        summary[i, "pvalue"] <= 0.01) {
                    count01 <- count01 + 1
                }
                if (summary[i, "median"] >= 0.2 & summary[i, "median"] < 0.3 &
                        summary[i, "pvalue"] <= 0.01) {
                    count02 <- count02 + 1
                }
                if (summary[i, "median"] >= 0.3 & summary[i, "median"] < 0.4 &
                        summary[i, "pvalue"] <= 0.01) {
                    count03 <- count03 + 1
                }
                if (summary[i, "median"] >= 0.4 & summary[i, "median"] < 0.5 &
                        summary[i, "pvalue"] <= 0.01) {
                    count04 <- count04 + 1
                }
                if (summary[i, "median"] >= 0.5 & summary[i, "median"] < 0.6 &
                        summary[i, "pvalue"] <= 0.01) {
                    count05 <- count05 + 1
                }
                if (summary[i, "median"] >= 0.6 & summary[i, "median"] < 0.7 &
                        summary[i, "pvalue"] <= 0.01) {
                    count06 <- count06 + 1
                }
                if (summary[i, "median"] >= 0.7 & summary[i, "median"] < 0.8 &
                        summary[i, "pvalue"] <= 0.01) {
                    count07 <- count07 + 1
                }
                if (summary[i, "median"] >= 0.8 & summary[i, "median"] < 0.9 &
                        summary[i, "pvalue"] <= 0.01) {
                    count08 <- count08 + 1
                }
                if (summary[i, "median"] >= 0.9 & summary[i, "median"] < 1 &
                        summary[i, "pvalue"] <= 0.01) {
                    count09 <- count09 + 1
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

makeBoxplot <- function(data, samples, score, bySamples=FALSE,
    col="standard", main="", xlab="Batch", ylab="Beta value", scoreCol=TRUE) {
    if (bySamples == FALSE) {
        ## get batch numbers
        batches <- sort(unique(samples$batch_id))
        ## prepare data for the box plot
        boxplotData <-list()
        for (i in 1:length(batches))
        {
            s<-as.character(samples[samples$batch_id == batches[i], 
                "sample_id"])
            boxplotData[[as.character(batches[i])]] <- as.vector(as.matrix(
                data[, s]))
            remove(s)
        }
        annotation <- rep("", length(batches))
        par(las=2)
        if (col == "standard") {
            coloring="yellow"
        }
        else {
            coloring=col
        }
        
        ## make boxplot
        boxplot(boxplotData, col=coloring, main=main, xlab=xlab, ylab=ylab,
            names=annotation, cex.axis=1, xaxt="n")

        if (scoreCol == TRUE) {
            ## get scoring colors for the x-axis labels
            colors <- c()
            for (i in batches) {
                if (score$BEscore[score$batch == i] >= 0.1) {
                    colors <- c(colors, "red")
                } else
                    if (score$BEscore[score$batch == i] >= 0.02 & 
                        score$BEscore[score$batch == i] < 1) {
                            colors <- c(colors, "orange")
                    } else
                        if (score$BEscore[score$batch == i] < 0.02) {
                            colors <- c(colors, "darkgreen")
                        }
            }
        }
        else {
            colors <- rep("black", length(batches))
        }

        ## add colored batch number to x-axis
        for (i in 1:length(batches)) {
            mtext(batches[i], side=1, line=1, at=i, las=2, cex=1,
                col=colors[i])
        }
    }
    else {
        ## get batch numbers
        batches <- sort(unique(samples$batch_id))
        ## boxplot separated by sample
        ## sort samples by batch
        samplesNew <- samples
        samplesNew <- samplesNew[order(samplesNew$batch_id), ]
        ## take sample_ids
        barcodes <- samplesNew$sample_id

        ## prepare data
        boxplotData<-list()
        for (i in 1:length(barcodes))
        {
            s <- as.character(samplesNew[samplesNew$sample_id == barcodes[i],
                "sample_id"])
            boxplotData[[as.character(barcodes[i])]] <- as.vector(as.matrix(
                data[, s]))
            remove(s)
        }

        annotation <- rep("", length(barcodes))
        par(las=2)
        ## get colors for the boxes
        batchesForColoring <- sort(samples$batch_id)
        ## store number of samples per batch
        countings <- as.vector(table(batchesForColoring))
        if (col == "standard") {
            coloring <- c()
            colnum <- 1
            for (i in 1:length(batches)) {
                if (colnum == 1) {
                    coloring <- c(coloring, rep("yellow", countings[i]))
                    colnum=2
                } else
                    if (colnum == 2) {
                        coloring <- c(coloring, rep("red", countings[i]))
                        colnum=3
                    } else
                        if (colnum == 3) {
                            coloring <- c(coloring, rep("blue", countings[i]))
                            colnum=4
                        } else
                            if (colnum == 4) {
                                coloring <- c(coloring, rep("green",
                                    countings[i]))
                                colnum=1
                            }
            }
        }
        else {
            coloring=col
        }

        ## make the boxplot
        boxplot(boxplotData, main=main, xlab=xlab, ylab=ylab,
            names=annotation, col=coloring, cex.axis=0.8, xaxt="n")

        ## get scoring colors for the x-axis labels
        if (scoreCol == TRUE) {
            colors <- c()
            for (i in batches) {
                if (score$BEscore[score$batch == i] >= 0.1) {
                    colors <- c(colors, "red")
                } else
                    if (score$BEscore[score$batch == i] >= 0.02 & 
                            score$BEscore[score$batch == i] < 1) {
                        colors <- c(colors, "orange")
                    } else
                        if (score$BEscore[score$batch == i] < 0.02) {
                            colors <- c(colors, "darkgreen")
                        }
            }
        }
        else {
            colors <- rep("black", length(batches))
        }

        ## calculate positions to add text for every batch
        positionsText <- c()
        positionCounter <- 0
        for (i in 1:length(batches)) {
            if (countings[i] == 1) {
                positionsText <- c(positionsText, positionCounter+1)
            } else
            if (countings[i] == 2) {
                positionsText <- c(positionsText, positionCounter+1)
            } else
            if (countings[i] == 3) {
                positionsText <- c(positionsText, positionCounter+2)
            } else
            ## nearly the mid position of the batch
            if (countings[i] >= 4) {
                positionsText <- c(positionsText, positionCounter+
                    ceiling(countings[i]/2))
            }
        positionCounter <- positionCounter + countings[i]
        }
        
        ## add colored batch number to x-axis
        batchCounter <- 1
        for (i in 1:positionCounter) {
            if (batchCounter == length(batches) + 1) {
                break;
            } else
            ## add number if i is on positionText of batch batchCounter
            if (i == positionsText[batchCounter]) {
                mtext(batches[batchCounter], side=1, line=1, at=i, las=2, 
                    cex=1, col=colors[batchCounter])
                batchCounter <- batchCounter + 1
            } else
            ## add "" if not
            if (i != positionsText[batchCounter]) {
                mtext("", side=1, line=1, at=i, las=2, cex=1)
            }
        }
    }
}

clearBEgenes <- function(data, samples, summary) {
    ## set beta values in data to NA for all found genes
    for (i in 1:nrow(summary))
    {
        data[summary$gene[i], as.character(
            samples$sample_id[samples$batch_id == summary$batch[i]])] <- NA
    }
    amount <- sum(is.na(data)) / (nrow(data) * ncol(data))
    message(paste(sum(is.na(data)), "values (", amount * 100,
        "% of the data) set to NA", sep=" "))
    return(data)
}

countValuesToPredict <- function(data) {
    amount <- sum(is.na(data)) / (nrow(data) * ncol(data))
    message(paste(sum(is.na(data)), "values (", amount * 100,
        "% of the data) set to NA", sep=" "))
    numberPredValues <- matrix(ncol=2, nrow=ncol(data))
    numberPredValues <- as.data.frame(numberPredValues)
    colnames(numberPredValues) <- c("sample", "num_pred_values")
    numberPredValues[, 1] <- colnames(data)
    numberPredValues[, 2] <- 0
    for (i in 1:nrow(data)) {
        for (j in 1:ncol(data)) {
            if (is.na(data[i, j])) {
                oldNumber <- numberPredValues[numberPredValues$sample == 
                    j, 2]
                newNumber <- oldNumber+1
                numberPredValues[numberPredValues$sample == j, 2] <- newNumber
            }
        }
    }
    return(numberPredValues=numberPredValues)
}

BEclear <- function(data, parallel=TRUE, cores=4, rowBlockSize=60,
    colBlockSize=60, epochs=50, outputFormat="RData", dir=getwd()) {

    ## needed functions

    ## returns vector with number of blocks, start - and stop position of
    ## last 2 blocks
    calcPositions <- function(num, blockSize) {
        rest <- num %% blockSize
        ## normal start and stop position if nrow can be divided by blockSize
        ## without a rest
        startPos1 <- num - (2*blockSize) - rest + 1
        startPos2 <- num - blockSize - rest + 1
        stopPos1 <- startPos2 - 1
        stopPos2 <- num
        ## special case if matrix has < blockSize rows
        if (num < blockSize) {
            startPos1 <- 0
            startPos2 <- 1
            stopPos1 <- 0
            blocks <- 1   
        } else {
            ## num can be divided by blockSize
            if ((num %% blockSize) == 0) {
                blocks <- num / blockSize # done
            } else {
                ## case not dividable by 100
                ## case rest < 75% of block size
                if (rest < floor(blockSize * 0.75)) {
                    blocks <- num %/% blockSize # done
                } else {
                    ## case rest > 75%
                    y <- blockSize + rest
                    blocks <- num %/% blockSize + 1 # done
                    ## case y dividable by 2
                    if (y %% 2 == 0) {
                        startPos1 <- num - y + 1
                        startPos2 <- startPos1 + (y %/% 2) + 1
                        stopPos1 <- startPos2 - 1    
                    }
                    else {
                        startPos1 <- num - y + 1
                        startPos2 <- startPos1 + (y %/% 2) + 1
                        stopPos1 <- startPos2 - 1
                    }
                }
            }
        }
        return <- c(blocks, startPos1, startPos2, stopPos1, stopPos2)
    }

    ## calculates every start and stop position for every block
    calcBlockFrame <- function(rowPositions, colPositions) {
        ## construct data frame with all row-and column star-and stop positions
        blocks <- data.frame(number=numeric(), rowStartPos=numeric(),
            rowStopPos=numeric(), colStartPos=numeric(),
            colStopPos=numeric())
        ## counters to choose the correct start - and stop positions within the
        ## loop
        rowCounter <- 1
        colCounter <- 1
        counter <- 1
        ## iteration over rows and columns to calculate the positions
        for (i in 0:(rowPos[1]-1)) {
            if (rowCounter < rowPos[1] - 1) {
                currentRowStart <- i*rowBlockSize+1
                currentRowStop <- i*rowBlockSize+rowBlockSize
                rowCounter <- rowCounter + 1
            } else {
                if (rowCounter == rowPos[1] - 1) {
                    currentRowStart <- rowPos[2]
                    currentRowStop <- rowPos[4]
                    rowCounter <- rowCounter + 1
                } else {
                    if (rowCounter == rowPos[1]) {
                        currentRowStart <- rowPos[3]
                        currentRowStop <- rowPos[5]
                    }
                }
            }
            for (j in 0:(colPos[1]-1)) {
                if (colCounter < colPos[1] - 1) {
                    currentColStart <- j*colBlockSize+1
                    currentColStop <- j*colBlockSize+colBlockSize
                    colCounter <- colCounter + 1
                } else {
                    if (colCounter == colPos[1] - 1) {
                        currentColStart <- colPos[2]
                        currentColStop <- colPos[4]
                        colCounter <- colCounter + 1
                    } else {
                        if (colCounter == colPos[1]) {
                            currentColStart <- colPos[3]
                            currentColStop <- colPos[5]
                            colCounter <- 1
                        }
                    }
                }
                blocks <- rbind(blocks, data.frame(number=counter,
                    rowStartPos=currentRowStart, rowStopPos=currentRowStop,
                    colStartPos=currentColStart, colStopPos=currentColStop))
                counter <- counter + 1
            }
        }
        return(blocks)
    }

    ## BEclear function
    BEclear <- function(block) {

        rowStartPosition <- blockFrame[[2]][[block]]
        rowStopPosition <- blockFrame[[3]][[block]]
        colStartPosition <- blockFrame[[4]][[block]]
        colStopPosition <- blockFrame[[5]][[block]]

        ## take one block of data (200 rows * 250 columns)
        D<-data[rowStartPosition:rowStopPosition,
            colStartPosition:colStopPosition]
        ## set NA data to 0
        D[is.na(D)] <- 0
        ## set D as matrix
        D <- as(as.matrix(D), "dgTMatrix")

        Dsummary <- summary(D)
        m <- nrow(D)               # number of rows
        n <- ncol(D)               # number of columns
        is <- Dsummary$i           # row of each revealed entry
        js <- Dsummary$j           # column of each revealed entry
        xs <- Dsummary$x           # value of each revealed entry
        N <- length(is)            # number of revealed entries

        ## counts for every index, the number of non-zero entries
        nnz <- function(indexes, l=max(indexes)) {
            ## array of l entries with 0
            n <- rep(0, l)
            ## get non-zero entries per index
            com <- tapply(indexes, indexes, length)
            ## save number of non-zero entries in n
            n[as.integer(names(com))] <- com
            return(n)
        }

        ## for each row, number of revealed entries=entries not zero
        nnzis <- nnz(is, m)
        ## for each column, number of revealed entries=entries not zero
        nnzjs <- nnz(js, n)
        ## changing values


        ## compute the loss of factorization LR
        loss <- function(L, R, lambda) {
            lossResult <- 0
            for (i in 1:N) {
                lossResult <- lossResult + ((D[is[i], js[i]] - 
                    (L[is[i], ]%*%R[, js[i]]))^2)
            }
            lossResult <- lossResult + lambda * (norm(L, type="F")^2)
            + lambda * (norm(R, type="F")^2)
            return(as.numeric(lossResult))
        }

        ## compute the local loss of the p-th revelead entry for
        ## factorization LR
        dlossp <- function(L, R, lambda, p) {
            i <- is[p]
            j <- js[p]
            x <- D[i, j]

            ## create two length-r vectors dLi and dRj
            ## * dLi[k] contains the gradient of the local loss with
            ## respect to L_ik
            ## * dRi[k] contains the gradient of the local loss with
            ## respect to R_kj
            dLi <- vector(mode="numeric", length=r)
            dRj <- vector(mode="numeric",  length=r)
            for(k in 1:r) {
                dLi[k] <- -2 * R[k, j] * (x - L[i, ]%*%R[, j]) + 2 * lambda * 
                    L[i, k]/nnzis[i]
                dRj[k] <- -2 * L[i, k] * (x - L[i, ]%*%R[, j]) + 2 * lambda * 
                    R[k, j]/nnzjs[j]
            }

            ## return the result as a list
            ## (elements can be accessed with x$Li or x$Ri, where x
            ## is the returned list)
            return(list(Li=dLi, Rj=dRj))
        }

        ## Run a gradient descent epoch (L and R are starting points,
        ## eps is stepsize)
        gdepoch <- function(L, R, lambda, eps) {
            ## create gradient matrices
            dL <- matrix(0, m, r)
            dR <- matrix(0, r, n)

            ## fill the gradient matrices using repeated calls to
            ## dlossp function
            for (i in 1:N) {
                dL[is[i], ] <- dL[is[i], ] + dlossp(L, R, lambda, i)$Li
                dR[, js[i]] <- dR[, js[i]] + dlossp(L, R, lambda, i)$Rj
            }

            ## perform a gradient step on L and R with step size eps
            ## by using the gradient matrices
            L <- L - eps * dL
            R <- R - eps * dR

            ## return result
            return(list(L=L, R=R))
        }

        ## Runner for gradient descent (or stochastic gradient descent) for the
        ## specified number of epoch
        runner <- function(f, L0r10, R0r10, lambda, epochs=epochs,
                eps=0.01) {
            LR <- list(L=L0r10, R=R0r10)
            curLoss <- loss(LR$L, LR$R, lambda)
            for (epoch in 1:epochs) {
                LR <- f(LR$L, LR$R, lambda, eps)

                ## bold driver step size update
                oldLoss <- curLoss
                curLoss <- loss(LR$L, LR$R, lambda)
                if (oldLoss < curLoss) { 
                    eps <- eps/2
                } else {
                    eps <- eps * 1.05
                }
            }
            return(LR)
        }


        r <- 10
        set.seed(1, kind="Mersenne-Twister")
        L0r10 <- matrix(rnorm(m * r), m, r) / sqrt(r)
        R0r10 <- matrix(rnorm(r * n), r, n) / sqrt(r)
        ## run LFM 
        resultGdr10l1 <- runner(gdepoch, L0r10, R0r10, 1, epochs=epochs,
            eps=0.01)
        
        dataTemp <- data[rowStartPosition:rowStopPosition,
            colStartPosition:colStopPosition]
        ## Predicted matrix
        D1 <- resultGdr10l1$L %*% resultGdr10l1$R
        ## Replace predicted values with known, keep predicted unknown values
        for(i1 in 1:nrow(dataTemp))
            for(j1 in 1:ncol(dataTemp))
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
            rowStopPosition, "col", colStartPosition, colStopPosition,
            "RData", sep=".")
        save(D1, file=paste(dir, filename, sep="/"))
        
        remove(D1, D, dataTemp)
        return(block)
    }

    ## parallel BEclear
    BEclearParallel <- function() {
        ## start parallel mode on 4 cpus
        sfInit(parallel=TRUE, cpus=cores, type="SOCK")
        ## load packages needed in parallel computing
        sfLibrary(Matrix)
        ## load variables needed in parallel computing
        sfExport("data", "blockFrame")
        ## load functions needed in parallel computing
        sfExport("BEclear")

        ## calculate BEclear in parallel mode
        result <- sfClusterApplyLB(blockNumbers, BEclear)

        ## stop parallel mode
        sfStop()

        return(unlist(result))
    }

    ## load and combine packages
    combineBlocks <- function(blockFrame, rowPos, colPos) {
        D1 <- NULL
        blockCounter <- 1
        for (i in 1:rowPos[1]) {
            rowStart <- blockFrame[[2]][[blockCounter]]
            rowStop <- blockFrame[[3]][[blockCounter]]
            for (j in 1:colPos[1]) {
                colStart <- blockFrame[[4]][[blockCounter]]
                colStop <- blockFrame[[5]][[blockCounter]]
                ## special case, only one column block
                if (colPos[1] == 1) {
                    blockName <- paste("D", blockCounter, sep="")
                    filename <- paste(blockName, "row", rowStart, rowStop, 
                        "col", colStart, colStop, "RData", sep=".")
                    file <- paste(dir, filename, sep="/")
                    load(file)
                    assign(paste(blockName, "row", rowStart, rowStop, "col",
                        colStart, colStop, sep="."), D1)
                    frame <- D1
                    currentRowGroup <- as.data.frame(frame)
                    blockCounter <- blockCounter + 1
                } else {
                    if (j == 1) {
                        blockName <- paste("D", blockCounter, sep="")
                        filename <- paste(blockName, "row", rowStart, rowStop,
                            "col", colStart, colStop, "RData", sep=".")
                        file <- paste(dir, filename, sep="/")
                        load(file)
                        assign(paste(blockName, "row", rowStart, rowStop,
                            "col", colStart, colStop, sep="."), D1)
                        frame <- D1
                        currentRowGroup <- as.data.frame(frame)
                        blockCounter <- blockCounter + 1
                    }
                    else {
                        blockName <- paste("D", blockCounter, sep="")
                        filename <- paste(blockName, "row", rowStart, rowStop,
                            "col", colStart, colStop, "RData", sep=".")
                        file <- paste(dir, filename, sep="/")
                        load(file)
                        assign(paste(blockName, "row", rowStart, rowStop,
                            "col", colStart, colStop, sep="."), D1)
                        frame <- D1
                        currentRowGroup <- cbind(currentRowGroup,
                            as.data.frame(frame))
                        blockCounter <- blockCounter + 1
                    }
                }
            }
            ## special case, only one row block
            if (rowPos[1] == 1) {
                predictedGenes <- currentRowGroup
            } else {
                if (i == 1) {
                    predictedGenes <- currentRowGroup
                }
                else {
                    predictedGenes <- rbind(predictedGenes, currentRowGroup)
                }
            }
        }
        return(as.matrix(predictedGenes))
    }

    ## start
    if (parallel == FALSE) {
        D1 <- NULL
        if (epochs<=0) {
            stop('number of epochs has to be greater than 0')
        } else if (cores<=0) {
            cores == 1
            stop('wrong number of cores, cores is set to one')
        } else if (rowBlockSize == 0 & colBlockSize == 0) {
            cat("BEclear is startet in non-parallel mode on the whole\n")
            cat("input matrix")
            ## start BEclear in non-parallel mode, just with one block
            blockFrame <- t(as.matrix(c(1, 1, nrow(data), 1, ncol(data))))
            blockFrame <- as.data.frame(blockFrame)
            blocksDone <- BEclear(1)
            ## load block
            blockName <- paste("D", blockFrame[1], sep="")
            filename <- paste(blockName, "row", blockFrame[2], blockFrame[3],
                "col", blockFrame[4], blockFrame[5], "RData", sep=".")
            file <- paste(dir, filename, sep="/")
            load(file)
            assign("predictedGenes", D1)
            
            ## remove stored single block file
            blockFilenames <- c(paste(blockName, "row", blockFrame[2],
                blockFrame[3], "col", blockFrame[4], blockFrame[5], "RData",
                sep="."))
            filedir <- paste(dir, blockFilenames, sep="/")
            file.remove(filedir)
            
            remove(D1, blockFrame)
            ## save block as predictedGenes
            predictedGenes <- as.data.frame(predictedGenes)
            if (outputFormat == "RData") {
                filename=paste("predicted.genes", "RData", sep=".")
                save(predictedGenes, file=paste(dir, filename, sep="/"))
            }
            if (outputFormat == "txt") {
                filename=paste("predicted.genes", "txt", sep=".")
                write.table(predictedGenes, 
                    file=paste(dir, filename, sep="/"),
                    row.names=TRUE, col.names=TRUE, sep="\t")
            }
            return(as.matrix(predictedGenes))    
        }
        else {
            cat("BEclear is startet in non-parallel mode;\n")
            cat(paste("block size: ", rowBlockSize, " x ", 
                colBlockSize, sep=""))
            ## run BEclear in non-parallel mode
            ## calculate start - and stop position for every block
            if (nrow(data) < rowBlockSize) {
                rowBlockSize <- nrow(data)
            }
            if (ncol(data) < colBlockSize) {
                colBlockSize <- ncol(data)
            }
            rowPos <- calcPositions(nrow(data), rowBlockSize)
            colPos <- calcPositions(ncol(data), colBlockSize)
            blockFrame <- calcBlockFrame(rowPos, colPos)
            remove(rowBlockSize, colBlockSize)

            ## only block numbers needed to distribute the blocks onto the
            ## worker slaves
            blockNumbers <- blockFrame[, 1]

            ## run BEclear in non-parallel mode
            for (i in blockNumbers) {
                blocksDone <- BEclear(i)
            }

            ## combine the blocks to the predictedGenes data.frame
            predictedGenes <- combineBlocks(blockFrame, rowPos, colPos)
            colnames(predictedGenes) <- colnames(data)
            
            ## remove all stored single block files
            blockFilenames <- c()
            for (i in 1:nrow(blockFrame)) {
                row <- paste("D", blockFrame$number[i], ".row.",
                    blockFrame$rowStartPos[i], ".", blockFrame$rowStopPos[i],
                    ".col.", blockFrame$colStartPos[i], ".",
                    blockFrame$colStopPos[i], ".RData", sep="")
                blockFilenames <- c(blockFilenames, row)
            }
            for (i in 1:length(blockFilenames)) {
                filedir <- paste(dir, blockFilenames[i], sep="/")
                file.remove(filedir)
            }
            
            remove(blockFrame, blockNumbers, colPos, rowPos)

            if (outputFormat == "RData") {                
                filename=paste("predicted.genes", "RData", sep=".")
                save(predictedGenes, file=paste(dir, filename, sep="/"))
            }
            if (outputFormat == "txt") {
                filename=paste("predicted.genes", "txt", sep=".")
                write.table(predictedGenes, 
                    file=paste(dir, filename, sep="/"),
                    row.names=TRUE, col.names=TRUE, sep="\t")
            }
            return(as.matrix(predictedGenes))
        }
    }
    else {
        D1 <- NULL
        if (epochs<= 0) {
            stop('number of epochs has to be greater than 0')
        } else {
            if (cores<=0) {
                cores == 1
                stop('wrong number of cores, cores is set to one')
            } else if (cores == 1) {
                if (rowBlockSize == 0 & colBlockSize == 0) {
                    cat("cores is set to one, BEclear is startet in \n")
                    cat("non-parallel mode on the whole input matrix")
                    ## start BEclear in non-parallel mode, just with one block
                    blockFrame <- t(as.matrix(c(1, 1, nrow(data), 1,
                        ncol(data))))
                    blockFrame <- as.data.frame(blockFrame)
                    blocksDone <- BEclear(1)
                    ## load block
                    blockName <- paste("D", blockFrame[1], sep="")
                    filename <- paste(blockName, "row", blockFrame[2],
                        blockFrame[3], "col", blockFrame[4], blockFrame[5],
                        "RData", sep=".")
                    file <- paste(dir, filename, sep="/")
                    load(file)
                    assign("predictedGenes", D1)
                    
                    ## remove stored single block file
                    blockFilenames <- c(paste(blockName, "row", blockFrame[2],
                        blockFrame[3], "col", blockFrame[4], blockFrame[5],
                        "RData", sep="."))
                    filedir <- paste(dir, blockFilenames, sep="/")
                    file.remove(filedir)
                    
                    remove(D1, blockFrame)
                    ## save block as predictedGenes
                    predictedGenes <- as.data.frame(predictedGenes)
                    if (outputFormat == "RData") {
                        filename=paste("predicted.genes", "RData", sep=".")
                        save(predictedGenes, 
                            file=paste(dir, filename, sep="/"))
                    }
                    if (outputFormat == "txt") {
                        filename=paste("predicted.genes", "txt", sep=".")
                        write.table(predictedGenes,
                            file=paste(dir, filename, sep="/"),
                            row.names=TRUE, col.names=TRUE, sep="\t")
                    }
                    return(as.matrix(predictedGenes))
                } else {
                    cat("cores is set to one, BEclear is startet in\n")
                    cat(paste("non-parallel mode; block size: ", rowBlockSize))
                    cat(paste(" x ", colBlockSize, sep=""))
                    ## run BEclear in non-parallel mode
                    ## calculate start - and stop position for every block
                    if (nrow(data) < rowBlockSize) {
                        rowBlockSize <- nrow(data)
                    }
                    if (ncol(data) < colBlockSize) {
                        colBlockSize <- ncol(data)
                    }
                    rowPos <- calcPositions(nrow(data), rowBlockSize)
                    colPos <- calcPositions(ncol(data), colBlockSize)
                    blockFrame <- calcBlockFrame(rowPos, colPos)
                    remove(rowBlockSize, colBlockSize)

                    ## only block numbers needed to distribute the blocks onto
                    ## the worker slaves
                    blockNumbers <- blockFrame[, 1]

                    ## run BEclear in non-parallel mode
                    for (i in blockNumbers) {
                        blocksDone <- BEclear(i)
                    }

                    ## combine the blocks to the predictedGenes data.frame
                    predictedGenes <- combineBlocks(blockFrame, rowPos,
                        colPos)
                    colnames(predictedGenes) <- colnames(data)
                    
                    ## remove all stored single block files
                    blockFilenames <- c()
                    for (i in 1:nrow(blockFrame)) {
                        row <- paste("D", blockFrame$number[i], ".row.",
                            blockFrame$rowStartPos[i], ".",
                            blockFrame$rowStopPos[i], ".col.",
                            blockFrame$colStartPos[i], ".",
                            blockFrame$colStopPos[i], ".RData", sep="")
                        blockFilenames <- c(blockFilenames, row)
                    }
                    for (i in 1:length(blockFilenames)) {
                        filedir <- paste(dir, blockFilenames[i], sep="/")
                        file.remove(filedir)
                    }
                    
                    remove(blockFrame, blockNumbers, colPos, rowPos)
                    if (outputFormat == "RData") {
                        filename=paste("predicted.genes", "RData", sep=".")
                        save(predictedGenes, 
                            file=paste(dir, filename, sep="/"))
                    }
                    if (outputFormat == "txt") {
                        filename=paste("predicted.genes", "txt", sep=".")
                        write.table(predictedGenes,
                            file=paste(dir, filename, sep="/"),
                            row.names=TRUE, col.names=TRUE, sep="\t")
                    }
                    return(as.matrix(predictedGenes))
                }
            } else {
                if (rowBlockSize == 0 & colBlockSize == 0) {
                    cat("block size is set to 0, BEclear is started in \n")
                    cat("non-parallel mode on the whole data matrix")
                    ## start BEclear in non-parallel mode, just with one block
                    blockFrame <- t(as.matrix(c(1, 1, nrow(data), 1,
                        ncol(data))))
                    blockFrame <- as.data.frame(blockFrame)
                    blocksDone <- BEclear(1)
                    ## load block
                    blockName <- paste("D", blockFrame[1], sep="")
                    filename <- paste(blockName, "row", blockFrame[2],
                        blockFrame[3], "col", blockFrame[4], blockFrame[5],
                        "RData", sep=".")
                    file <- paste(dir, filename, sep="/")
                    load(file)
                    assign("predictedGenes", D1)
                    
                    ## remove stored single block file
                    blockFilenames <- c(paste(blockName, "row", blockFrame[2],
                        blockFrame[3], "col", blockFrame[4], blockFrame[5],
                        "RData", sep="."))
                    filedir <- paste(dir, blockFilenames, sep="/")
                    file.remove(filedir)
                    
                    remove(D1, blockFrame)
                    ## save block as predictedGenes
                    predictedGenes <- as.data.frame(predictedGenes)
                    if (outputFormat == "RData") {
                        filename=paste("predicted.genes", "RData", sep=".")
                        save(predictedGenes, 
                            file=paste(dir, filename, sep="/"))
                    }
                    if (outputFormat == "txt") {
                        filename=paste("predicted.genes", "txt", sep=".")
                        write.table(predictedGenes,
                            file=paste(dir, filename, sep="/"),
                            row.names=TRUE, col.names=TRUE, sep="\t")
                    }
                    return(as.matrix(predictedGenes))
                } else {
                    ## run BEclear in parallel mode
                    cat("BEclear is startet in parallel mode; block size:\n")
                    cat(paste(rowBlockSize, " x ", colBlockSize, sep=""))
                    ## calculate start - and stop position for every block
                    if (nrow(data) < rowBlockSize) {
                        rowBlockSize <- nrow(data)
                    }
                    if (ncol(data) < colBlockSize) {
                        colBlockSize <- ncol(data)
                    }
                    rowPos <- calcPositions(nrow(data), rowBlockSize)
                    colPos <- calcPositions(ncol(data), colBlockSize)
                    blockFrame <- calcBlockFrame(rowPos, colPos)
                    remove(rowBlockSize, colBlockSize)

                    ## only block numbers needed to distribute the blocks onto 
                    ## the worker slaves
                    blockNumbers <- blockFrame[, 1]
                    ## run BEclear in parallel mode
                    blocksDone <- BEclearParallel()

                    ## combine the blocks to the predictedGenes data.frame
                    predictedGenes <- combineBlocks(blockFrame, rowPos,
                        colPos)
                    colnames(predictedGenes) <- colnames(data)
                    
                    ## remove all stored single block files
                    blockFilenames <- c()
                    for (i in 1:nrow(blockFrame)) {
                        row <- paste("D", blockFrame$number[i], ".row.",
                            blockFrame$rowStartPos[i], ".",
                            blockFrame$rowStopPos[i], ".col.",
                            blockFrame$colStartPos[i], ".",
                            blockFrame$colStopPos[i], ".RData", sep="")
                        blockFilenames <- c(blockFilenames, row)
                    }
                    for (i in 1:length(blockFilenames)) {
                        filedir <- paste(dir, blockFilenames[i], sep="/")
                        file.remove(filedir)
                    }
                    
                    remove(blockFrame, blockNumbers, colPos, rowPos)
                    if (outputFormat == "RData") {
                        filename=paste("predicted.genes", "RData", sep=".")
                        save(predictedGenes, 
                            file=paste(dir, filename, sep="/"))
                    }
                    if (outputFormat == "txt") {
                        filename=paste("predicted.genes", "txt", sep=".")
                        write.table(predictedGenes,
                            file=paste(dir, filename, sep="/"),
                            row.names=TRUE, col.names=TRUE, sep="\t")
                    }
                    return(as.matrix(predictedGenes))
                }
            }
        }
    }
}

findWrongValues <- function(data) {
    ## find entries > 1 or < 0 in the matrix, store wrong entries in a
    ## data.frame
    wrongEntries <- as.data.frame(matrix(nrow=0, ncol=4))
    colnames(wrongEntries) <- c("level", "row", "col", "value")
    counter <- 1
    for(i in 1:nrow(data)) {
        for(j in 1:ncol(data)) {
            if(data[i, j] > 1) {
                wrongEntries[counter, "level"] <- 1
                wrongEntries[counter, "row"] <- i
                wrongEntries[counter, "col"] <- j
                wrongEntries[counter, "value"] <- data[i, j]
                counter <- counter + 1
            }
            if(data[i, j] < 0) {
                wrongEntries[counter, "level"] <- 0
                wrongEntries[counter, "row"] <- i
                wrongEntries[counter, "col"] <- j
                wrongEntries[counter, "value"] <- data[i, j]
                counter <- counter + 1
            }
        }
    }
    cat(paste(counter-1, " wrong values found", sep=""))
    return(wrongEntries)
}

replaceWrongValues <- function(data) {
    counter <- 0
    for(i in 1:nrow(data)) {
        for(j in 1:ncol(data)) {
            if(data[i, j] > 1) {
                data[i, j] <- 1
                counter <- counter + 1
            }
            if(data[i, j] < 0) {
                data[i, j] <- 0
                counter <- counter + 1
            }
        }
    }
    cat(paste(counter, " wrong values replaced", sep=""))
    return(data)
}


correctBatchEffect <- function(data, samples, parallel=TRUE, cores=4,
    adjusted=TRUE, method="fdr", rowBlockSize=60, colBlockSize=60, epochs=50,
    outputFormat="RData", dir=getwd()) {

    med <- calcMedians(data, samples, parallel, cores)
    pval <- calcPvalues(data, samples, parallel, cores, adjusted, method)
    sum <- calcSummary(med, pval)
    score <- calcScore(data, samples, sum)
    cleared <- clearBEgenes(data, samples, sum)
    predicted <- BEclear(data, parallel, cores, rowBlockSize, colBlockSize,
        epochs, outputFormat, dir)
    corrected <- replaceWrongValues(predicted)

    return(list(medians=med, pvals=pval, summary=sum, scoreTable=score,
        clearedData=cleared, predictedData=predicted,
        correctedPredictedData=corrected))
}
