\name{data_mix_stratum}
\alias{data_mix_stratum}
\docType{data}
\title{
An example with a mixture of endpoint types with three strata.
}
\description{
This is a simulated data example with a mixture of two continuous and one time-to-event endpoints with three strata. The first endpoint is time-to-event and the second and third endpoints are continuous.
}
\usage{data("data_mix_stratum")}
\format{
  A data frame with 400 observations on the following 6 variables.
  \describe{
    \item{\code{arm}}{A vector for the treatment groups.}
    \item{\code{stratum}}{A vector for the stratum.}
    \item{\code{Delta_1}}{A vector for the event status of the first endpoint (1=event, 0=censored).}
    \item{\code{Y_1}}{A vector for the outcome of the first endpoint.}
    \item{\code{Y_2}}{A vector for the outcome of the second endpoint.}
    \item{\code{Y_3}}{A vector for the outcome of the third endpoint.}
  }
}
\examples{
data(data_mix_stratum)
str(data_mix_stratum)
}
\keyword{datasets}
