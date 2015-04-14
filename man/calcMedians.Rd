\name{calcMedians}
\alias{calcMedians}
\title{Calculate median difference values}
\description{
Compares the median value of all beta values belonging to one batch with the
median value of all beta values belonging to all other batches. Returns a
matrix containing this median difference value for every gene in every batch,
columns define the batch numbers, rows the gene names.
}
\usage{
calcMedians(data, samples, parallel=TRUE, cores=4)
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
}

\value{
a matrix containing median comparison values for all genes in all batches, the
column names define the batch numbers, row names are the same gene names as
contained in the input matrix.
}

\seealso{
\code{\link{snowfall}}
\code{\link{correctBatchEffect}}
}
\examples{

## Shortly running example. For a more realistic example that takes
## some more time, run the same procedure with the full BEclearData
## dataset.

## Calculate median comparison values in non-parallel mode
data(BEclearData)
ex.data <- ex.data[31:90,7:26]
ex.samples <- ex.samples[7:26,]
medians <- calcMedians(data=ex.data, samples=ex.samples, parallel=FALSE)

}

\keyword{ median comparison }
