\name{calcPvalues}
\alias{calcPvalues}
\title{Calculate p-values}
\description{
Compares the distribution of all beta values corresponding to one batch with
the distribution of all beta values corresponding to all other batches and
retuns a p-value which defines if the distributions are the same or not.
Standard two sided Kolmogorov-Smirnov test is used to calculate the (adjusted)
p-values.
}
\usage{
calcPvalues(data, samples, parallel=TRUE, cores=4, adjusted=TRUE, method="fdr")
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
    function in parallel mode.}
    \item{adjusted}{should the p-values be adjusted or not, see "method" for
    available adjustment methods.}
    \item{method}{adjustment method for p-value adjustment (if TRUE), default
    method is "false discovery rate adjustment", for other available methods
    see the description of the used standard R package \code{\link{p.adjust}}.}
}

\value{
a matrix containing p-values for all genes in all batches, the column names
define the batch numbers, row names are the same gene names as contained in the
input matrix.

}

\seealso{
\code{\link{snowfall}}
\code{\link{ks.test}}
\code{\link{p.adjust}}
\code{\link{correctBatchEffect}}
}
\examples{
## Shortly running example. For a more realistic example that takes
## some more time, run the same procedure with the full BEclearData
## dataset.

## Calculate fdr-adjusted p-values in non-parallel mode
data(BEclearData)
ex.data <- ex.data[31:90,7:26]
ex.samples <- ex.samples[7:26,]
pvals <- calcPvalues(data=ex.data, samples=ex.samples, parallel=FALSE,
    adjusted=TRUE, method="fdr")
}

\keyword{ p-value }
\keyword{ p-value adjustment }