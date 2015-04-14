\name{BEclear example data}
\alias{BEclearData}
\alias{BEclear example data}
\alias{ex.data}
\alias{ex.samples}
\docType{data}
\title{Example data set for the BEclear-package}
\description{
Example data needed to run all examples contained in the documentation files of
the BEclear-package
}
\usage{data(BEclearData)}
\format{
Contains a data matrix and a samples data frame which can be used to run the
examples of the BEclear package.
    \describe{
        \item{\code{ex.data}}{An example data matrix that is filled with
        beta values originally stemming from breast cancer data from the
        TCGA portal [1], colnames are sample ids, rownames are gene names.
        Generally, beta values are claculated by dividing the methylated signal
        by the sum of the unmethylated and methylated signals from a DNA
        methylation microrarray. The sample data used here contains averaged
        beta values of probes that belong to promoter regions of single genes.
        Another possibility would be to use beta values of single probes,
        whereby the probe names should then be used instead of the gene names
        as rownames of the matrix.
        }
        \item{\code{ex.samples}}{An example data frame containing a column
        for the sample id and a column for the corresponding batch id.}
    }
}
\value{
    A data frame (ex.samples) and a matrix (ex.data)
}
\references{
[1] http://cancergenome.nih.gov/
}
\examples{
\dontrun{
data(BEclearData)

## Whole procedure that has to be done to make a prediction using the BEclear
## function

# Calculate median difference values and p-values
meds <- calcMedians(data=ex.data, samples=ex.samples)
pvals <- calcPvalues(data=ex.data, samples=ex.samples)

# Summarize p-values and median differences for batch affected genes
sum <- calcSummary(medians=meds, pvalues=pvals)

# Set entries defined by the summary to NA
clearedMatrix <- clearBEgenes(data=ex.data, samples=ex.samples, summary=sum)

# Predict the missing entries with standard parameters
predicted <- BEclear(data=clearedMatrix) 
}
}
\keyword{datasets}
\keyword{data}
\keyword{BEclear}