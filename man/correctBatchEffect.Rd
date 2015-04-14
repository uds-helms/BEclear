\name{correctBatchEffect}
\alias{correctBatchEffect}
\title{Correct a batch effect in DNA methylation data}
\description{
This method combines most functions of the \code{\link{BEclear-package}} to
one. The method performs the whole process of searching for batch effects and
automatically correct them for a matrix of beta values stemming from DNA
methylation data.
}
\usage{
correctBatchEffect(data, samples, parallel=TRUE, cores=4, adjusted=TRUE,
    method="fdr", rowBlockSize=60, colBlockSize=60, epochs=50,
    outputFormat="RData", dir=getwd())
}
\arguments{
    \item{data}{any matrix filled with beta values, column names have to be
    sample_ids corresponding to the ids listed in "samples", row names have to
    be gene names.}
    \item{samples}{data frame with two columns, the first column has to contain
    the sample numbers, the second column has to contain the corresponding
    batch number. Colnames have to be named as "sample_id" and "batch_id".}
    \item{parallel}{should the calculation be done in parallel mode? 
    \code{\link{snowfall}} package is needed to run the function in parallel
    mode.}
    \item{cores}{if running in parallel mode, define the number of cores used
    for the calculation. \code{\link{snowfall}} package is needed to run the
    function in parallel mode}
    \item{adjusted}{should the p-values be adjusted or not, see "method" for
    available adjustment methods.}
    \item{method}{adjustment method for p-value adjustment (if TRUE), default
    method is "false discovery rate adjustment", for other available methods
    see the description of the used standard R package \code{\link{p.adjust}}.
    See \code{\link{calcPvalues}} for more information}   
    \item{rowBlockSize}{the number of rows that is used in a block if the
    function is run in parallel mode and/or not on the whole matrix. Set this,
    and the "colBlockSize" parameter to 0 if you want to run the function on
    the whole input matrix. See \code{\link{BEclear}} and especially the
    details section of See \code{\link{BEclear}} for more information about
    this feature.}
    \item{colBlockSize}{the number of columns that is used in a block if the
    function is run in parallel mode and/or not on the whole matrix. Set this,
    and the "rowBlockSize" parameter to 0 if you want to run the function on
    the whole input matrix. See \code{\link{BEclear}} and especially the
    details section of See \code{\link{BEclear}} for more information about
    this feature.}
    \item{epochs}{the number of iterations used in the gradient descent
    algorithm to predict the missing entries in the data matrix. See
    \code{\link{BEclear}} for more information.}
    \item{outputFormat}{you can choose if the finally returned data matrix
    should be saved as an .RData file or as a tab-delimited .txt file in the
    specified directory. Allowed values are "RData" and "txt". See
    \code{\link{BEclear}} for more information.}
    \item{dir}{set the path to a directory the predicted matrix should be
    stored. The current working directory is defined as default parameter.}
}
\details{
The function performs the whole process of searching for batch effects and
automatically correct them for a matrix of beta values stemming from DNA
methylation data. Thereby, the function is running most of the functions from
the \code{\link{BEclear-package}} in a logical order.

First, median comparison values are calculated by the \code{\link{calcMedians}}
function, followed by the calculation of p-values by the
\code{\link{calcPvalues}} function.
With the results from the median comparison and p-value calculation, a summary
data frame is build using the \code{\link{calcSummary}} function, and a scoring
table is established by the \code{\link{calcScore}} function.
Now, found entries from the summary are set to NA in the input matrix using the
\code{\link{clearBEgenes}} function, then the \code{\link{BEclear}} function is
used to predict the missing values and at the end, possibly existing wrongly
predicted entries (values lower than 0 or greater than 1) are corrected using
the \code{\link{replaceWrongValues}} function.
}
\value{
A list containing the following fields (for detailed information look at the
documentations of the corresponding functions):
\item{medians}{A data.frame containing all median comparison values
corresponding to the input matrix.}
\item{pvalues}{A data.frame containing all p-values corresponding to the input
matrix.}
\item{summary}{The summarized results of the median comparison and p-value
calculation.}
\item{score.table}{A data.frame containing the number of found genes and a
BEscore for every batch.}
\item{cleared.data}{the input matrix with all values defined in the summary
set to NA.}
\item{predicted.data}{the input matrix after all previously NA values have been
predicted.}
\item{corrected.predicted.data}{the predicted matrix after the correction for
wrongly predicted values.}
}
\references{
Y. Koren, R. Bell, C. Volinsky, Matrix factorization techniques for recommender
systems, IEEE Computer, 42(8), p. 30-37, 2009,
http://research.yahoo.com/pub/2859

E. Candes, B. Recht, Exact matrix completion via convex optimization,
Communications of the ACM, 55(6), p. 111-119, 2012,
http://doi.acm.org/10.114/2184319.2184343
}

\seealso{
\code{\link{calcMedians}}
\code{\link{calcPvalues}}
\code{\link{calcSummary}}
\code{\link{calcScore}}
\code{\link{clearBEgenes}}
\code{\link{BEclear}}
\code{\link{replaceWrongValues}}
}
\examples{

## Shortly running example. For a more realistic example that takes
## some more time, run the same procedure with the full BEclearData
## dataset.

## Whole procedure that has to be done to use this function.
## Correct the example data for a batch effect
data(BEclearData)
ex.data <- ex.data[31:90,7:26]
ex.samples <- ex.samples[7:26,]

# Note that row- and block sizes are just set to 10 to get a short runtime. To
# use these parameters, either use the default values or please note the
# description in the details section of \code{\link{BEclear}}
result <- correctBatchEffect(data=ex.data, samples=ex.samples, parallel=TRUE,
    cores=4, adjusted=TRUE, method="fdr", rowBlockSize=10, colBlockSize=10,
    epochs=50, outputFormat="RData", dir=getwd())

# Unlist variables
medians <- result$medians
pvals <- result$pvals
summary <- result$summary
score <- result$score.table
cleared <- result$clearedData
predicted <- result$predictedData
corrected <- result$correctedPredictedData

}

\keyword{ BEclear }
\keyword{ batch effect correction }
\keyword{ DNA methylation }
\keyword{ median comparison }
\keyword{ p-value }
