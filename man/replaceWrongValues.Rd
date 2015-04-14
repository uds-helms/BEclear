\name{replaceWrongValues}
\alias{replaceWrongValues}
\title{Replace DNA methylation values out of the boundaries}
\description{A method which replaces values below 0 or beyond 1 contained in
the input matrix. These wrong entries are replaced by 0 or 1, respectively.
Note that this method is especially designed for DNA methylation data.}
\usage{
replaceWrongValues(data)
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
replaces these values either by 0 or 1, respectively. Another method called
\code{\link{findWrongValues}} returns a list of existing wrong values and can
be run before the replacement.
}
\value{
Returns the input matrix with every value previously below 0 changed to 0 and
every value previously beyond 1 changed to 1.
}

\seealso{
\code{\link{BEclear}}
\code{\link{findWrongValues}}
\code{\link{correctBatchEffect}}
}
\examples{
data(BEclearCorrected)
# Replace wrongly predicted values
corrected <- replaceWrongValues(data=ex.corrected.data)

}
\keyword{ BEclear }
\keyword{ batch effect correction }
\keyword{ wrong values }