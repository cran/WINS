\name{data_mix}
\alias{data_mix}
\docType{data}
\title{
An example with a mixture of endpoint types.
}
\description{
This is a simulated data example with a mixture of two continuous and one time-to-event endpoints. The first endpoint is time-to-event and the second and third endpoints are continuous.
}
\usage{data("data_mix")}
\format{
  A data frame with 400 observations on the following 5 variables.
  \describe{
    \item{\code{id}}{A vector for the patient id.}
    \item{\code{arm}}{A vector for the treatment groups.}
    \item{\code{Delta_1}}{A vector for the event status of the first endpoint (1=event, 0=censored).}
    \item{\code{Y_1}}{A vector for the outcome of the first endpoint.}
    \item{\code{Y_2}}{A vector for the outcome of the second endpoint.}
    \item{\code{Y_3}}{A vector for the outcome of the third endpoint.}
  }
}
\examples{
data(data_mix)
str(data_mix)
}
\keyword{datasets}
