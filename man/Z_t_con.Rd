\name{Z_t_con}
\alias{Z_t_con}
\docType{data}
\title{
Covariate history in the control group.
}
\description{
This is a simulated data example for the covariate history in the control group.
}
\usage{data("data_tte")}
\format{
  A data frame with 796 observations on the following 4 variables.
  \describe{
    \item{\code{id}}{A vector for the patient id.}
    \item{\code{time}}{A vector for the observed time of the covariate.}
    \item{\code{Z1}}{A vector for the observed value of the first covariate.}
    \item{\code{Z2}}{A vector for the observed value of the second covariate.}
  }
}
\examples{
data(data_tte)
str(Z_t_con)
}
\keyword{datasets}
