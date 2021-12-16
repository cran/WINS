\name{data_tte}
\alias{data_tte}
\docType{data}
\title{
An example with three TTE endpoints.
}
\description{
This is a simulated data example with three time-to-event endpoints.
}
\usage{data("data_tte")}
\format{
  A data frame with 400 observations on the following 7 variables.
  \describe{
    \item{\code{arm}}{A vector for the treatment groups.}
    \item{\code{Delta_1}}{A vector for the event status of the first endpoint (1=event, 0=censored).}
    \item{\code{Delta_2}}{A vector for the event status of the second endpoint (1=event, 0=censored).}
    \item{\code{Delta_3}}{A vector for the event status of the third endpoint (1=event, 0=censored).}
    \item{\code{Y_1}}{A vector for the outcome of the first endpoint.}
    \item{\code{Y_2}}{A vector for the outcome of the second endpoint.}
    \item{\code{Y_3}}{A vector for the outcome of the third endpoint.}
  }
}
\examples{
data(data_tte)
str(data_tte)
}
\keyword{datasets}
