\name{data_continuous}
\alias{data_continuous}
\docType{data}
\title{
An example with three continuous endpoints.
}
\description{
This is a simulated data example with three continuous endpoints.
}
\usage{data("data_continuous")}
\format{
  A data frame with 250 observations on the following 4 variables.
  \describe{
    \item{\code{arm}}{A vector for the treatment groups.}
    \item{\code{Y_1}}{A vector for the outcome of the first endpoint.}
    \item{\code{Y_2}}{A vector for the outcome of the second endpoint.}
    \item{\code{Y_3}}{A vector for the outcome of the third endpoint.}
  }
}
\examples{
data(data_continuous)
str(data_continuous)
}
\keyword{datasets}
