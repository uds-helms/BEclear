\name{BEclear-package}
\alias{BEclear-package}
\alias{BEclearCorrected}
\docType{package}
\title{
Correct for batch effects in DNA methylation data
}
\description{
Provides some functions to detect and correct for batch effects in DNA
methylation data.
The core function \code{\link{BEclear}} is based on Latent Factor Models and
can also be used to predict missing values in any other matrix containing real
numbers.
}
\details{
\tabular{ll}{
Package: \tab BEclear\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-03\cr
License: \tab GPL-2\cr
}
\code{\link{correctBatchEffect}}:
The function combines most functions of the \code{\link{BEclear-package}} to
one. This function performs the whole process of searching for batch effects
and automatically correct them for a matrix of beta values stemming from DNA
methylation data.

\code{\link{BEclear}}:
This function predicts the missing entries of an input matrix (NA values)
through the use of a Latent Factor Model.

\code{\link{calcMedians}}:
Compares the median value of all beta values belonging to one batch with the
median value of all beta values belonging to all other batches. Returns a 
matrix containing this median difference value for every gene in every batch,
columns define the batch numbers, rows the gene names.

\code{\link{calcPvalues}}:
Compares the distribution of all beta values corresponding to one batch with
the distribution of all beta values corresponding to all other batches and
retuns a p-value which defines if the distributions are the same or not.

\code{\link{calcSummary}}:
Summarizes the results of the median comparison function
\code{\link{calcMedians}} and the p-value calculation function
\code{\link{calcPvalues}}. Should be used with the matrices originating from
these two functions.

\code{\link{calcScore}}:
Returns a table with the number of found genes with found p-values less or
equal to 0.01 and median values greater or equal to 0.05. A score is calculated
depending on the number of found genes as well as the magnitude of the median
difference values, this score is divided by the overall number of genes in the
data and returned as "BEscore". See the methods details for further information
and details about the score calculation.

\code{\link{makeBoxplot}}:
A simple \code{\link{boxplot}} is done with boxes either separated by batches
or by samples and describe the five number summary of all beta values
corresponding to a batch or a sample, respectively. The batch_ids are shown on
the x-axis with a coloring corresponding to the BEscore.

\code{\link{clearBEgenes}}:
A function that simply sets all values to NA which were previously found by
median value comparison and p-value calculation and are stored in a summary.
The summary defines which values in the data matrix are set to NA.

\code{\link{countValuesToPredict}}:
Simple function that counts all values in a matrix which are NA

\code{\link{findWrongValues}}:
A method which lists values below 0 or beyond 1 contained in the input matrix.
The wrong entries are stored in a data.frame together with the corresponding
row and column position of the matrix.

\code{\link{replaceWrongValues}}:
A method which replaces values below 0 or beyond 1 contained in the input
matrix. These wrong entries are replaced by 0 or 1, respectively.

}
\author{
Ruslan Akulenko, Markus Merl

Maintainer: Markus Merl <merl.markus@googlemail.com>
}
\references{
Y. Koren, R. Bell, C. Volinsky, Matrix factorization techniques for recommender
systems, IEEE Computer, 42(8), p. 30-37, 2009,
http://research.yahoo.com/pub/2859

E. Candes, B. Recht, Exact matrix completion via convex optimization,
Communications of the ACM, 55(6), p. 111-119, 2012, 
http://doi.acm.org/10.114/2184319.2184343
}

\keyword{ BEclear }
\keyword{ batch effect correction }
\keyword{ DNA methylation }
\examples{
data(BEclearData)

## Calculate median comparison values in non-parallel mode
med <- calcMedians(data=ex.data, samples=ex.samples, parallel=FALSE)

## Calculate fdr-adjusted p-values in non-parallel mode
pvals <- calcPvalues(data=ex.data, samples=ex.samples, parallel=FALSE,
    adjusted=TRUE, method="fdr")

## Summarize p-values and median differences for batch affected genes
sum <- calcSummary(medians=med, pvalues=pvals)

## Calculates the score table
score.table <- calcScore(data=ex.data, samples=ex.samples, summary=sum)

## Simple boxplot for the example data separated by batch
makeBoxplot(data=ex.data, samples=ex.samples, score=score.table,
    bySamples=FALSE, main="Some box plot")

## Simple boxplot for the example data separated by samples
makeBoxplot(data=ex.data, samples=ex.samples, score=score.table,
    bySamples=TRUE, main="Some box plot")

## Sets assumed batch affected entries to NA
cleared <- clearBEgenes(data=ex.data, samples=ex.samples, summary=sum)

## Counts and stores number of entries to predict
numberOfEntries <- countValuesToPredict(data=cleared)

\dontrun{
## Predicts the missing entries
predicted <- BEclear(data=cleared)

## Find wrongly predicted entries
wrongEntries <- findWrongValues(data=predicted)

## Replace wrongly predicted entries
corrected <- replaceWrongValues(data=predicted)
}
}
