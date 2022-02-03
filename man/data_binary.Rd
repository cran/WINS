\name{data_binary}
\alias{data_binary}
\docType{data}
\title{
An example with three binary endpoints.
}
\description{
This is a simulated data example with three binary endpoints.
}
\usage{data("data_binary")}
\format{
  A data frame with 250 observations on the following 4 variables.
  \describe{
    \item{\code{id}}{A vector for the patient id.}
    \item{\code{arm}}{A vector for the treatment groups.}
    \item{\code{Y_1}}{A vector for the outcome of the first endpoint.}
    \item{\code{Y_2}}{A vector for the outcome of the second endpoint.}
    \item{\code{Y_3}}{A vector for the outcome of the third endpoint.}
  }
}
\examples{
data(data_binary)
str(data_binary)
}
\keyword{datasets}
