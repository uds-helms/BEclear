\name{findWrongValues}
\alias{findWrongValues}
\title{Find DNA methylation values out of the boundaries}
\description{A method which lists values below 0 or beyond 1 contained in the
input matrix. The wrong entries are stored in a data.frame together with the
corresponding row and column position of the matrix. Note that this method is
especially designed for DNA methylation data.}
\usage{
findWrongValues(data)
}
\arguments{
    \item{data}{any matrix filled with values that normally should be bounded
    between 0 and 1.}
}
\details{
Note that this method is especially designed to run after the batch effect
correction of DNA methylation data, e.g. with the \code{\link{BEclear}} method.
It can happen, that the predicted values are lying slightly below the lower
bound of 0 or beyond the upper bound of 1. This method finds these inaccurately
predicted entries. Another method called \code{\link{replaceWrongValues}}
replaces these values either by 0 or 1, respectively.
}
\value{
A data frame containing the columns "level", "row", "col" and "value" defining
if the wrong value is below 0 or beyond 1 (level = 0 or level = 1), the row
position and the column position in the input matrix and the wrong value
itself, respectively.
}

\seealso{
\code{\link{BEclear}}
\code{\link{replaceWrongValues}}
\code{\link{correctBatchEffect}}
}
\examples{
data(BEclearCorrected)
# Find wrongly predicted values
wrongEntries <- findWrongValues(data=ex.corrected.data)

}
\keyword{ BEclear }
\keyword{ batch effect correction }
\keyword{ wrong values }