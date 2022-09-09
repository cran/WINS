\name{partition_t.plot}
\alias{partition_t.plot}
\title{
Plot the Win Proportion over the Study Time.
}
\description{
A plot for the win proportions of the treatment/control group over the study time, for time-to-event endpoints only. The function "win.stat" is called to calculate the win proportions at each time in the plot function.
}
\usage{
partition_t.plot(data, Ctime = Inf, arm.name = c(1, 2), priority = c(1, 2),
censoring_adjust = "No", Z_t_trt = NULL, Z_t_con = NULL, tau = 0,
plotTimeUnit = NULL, trt_group = c("both","trt","con"),
win.strategy = NULL, ...)
}
\arguments{
  \item{data}{The analysis dataset which contains the following variables:
    \itemize{
      \item{arm: }{A vector for the treatment groups. Alternative names for "arm" include "trt", "treat" and "treatment".}
      \item{id: }{A vector for the patient id.}
      \item{stratum: }{A vector for the stratum.}
      \item{Delta_j: }{A vector for the event status of the j-th endpoint if the endpoint is a time-to-event outcome (1=event, 0=censored).}
      \item{Y_j: }{A vector for the outcome of the j-th endpoint, for time-to-event outcome, Y_j is a vector for the observed time.}
      \item{Start_time: } A vector for the time when each of the individuals is first accrued to study.
    }
  }
  \item{Ctime}{A vector of study times, at which the win proportions are to be plotted. Study time for a subject is the time since the subject is accrued, default is Inf.}
  \item{arm.name}{A vector for the two treatments, default to be c(1,2). The first value is for the treatment group, and the second value is for the control group.}
  \item{priority}{Importance order (from the most to the least important). For example, given three endpoints with the importance order as Endpoint 3, Endpoint 2, and Endpoint 1, input priority = c(3,2,1).}
  \item{censoring_adjust}{The method to adjust censoring for the kernal functions. Possible choices are listed below.
    \itemize{
      \item{"No": } {Without an adjustment.}
      \item{"IPCW": } {IPCW (inverse-probability-of-censoring weighting) adjustment for censoring using the Kaplan-Meier estimator.}
      \item{"CovIPCW": } {IPCW adjustment for censoring using the time dependent Cox model.}
    }
  }
  \item{Z_t_trt}{A matrix of the covariate history in the treatment group, each row is a (p+2) vector: the first two columns are id, time, the other p columns are the covariates (vector of length p). The baseline covariates are provided as the observed covariates corresponding to time 0.}
  \item{Z_t_con}{A matrix of the covariate history in the control group, each row is a (p+2) vector: the first two columns are id, time, the other p columns are the covariates (vector of length p). The baseline covariates are provided as the observed covariates corresponding to time 0.}
  \item{tau}{A vector of numerical value for the magnitude of difference to determine win/loss/tie for each endpoint. Tau is applicable for TTE endpoints and continuous endpoints; tau is fixed as 0 for binary endpoints. Default is 0 for all endpoints.}
  \item{plotTimeUnit}{The time units, e.g., "days", "months", "years", default as NULL.}
  \item{trt_group}{An argument to decide the arms for which the win proportions over time are plotted, possible choices include "both", "trt" and "con". Default is "both", i.e., by default the win proportion over time is plotted both for the treatment group and for the control group.}
  \item{win.strategy}{The strategy to determine the win status. Default as NULL. If NULL, the default win strategy funtion "win.strategy.default" is called, see \link{win.strategy.default} for more details. Users can also define their own "win.strategy" function.}
  \item{...}{Argument passed from user defined functions "win.strategy" if there is any. For instructions on this "win.strategy" function, see \link{win.strategy.default} for more details.}
}
\value{
A ggplot2 object.
}
\examples{
\donttest{
#### An simulated example with two TTE endpoints.
data <- sim.data(n_trt = 200, n_con = 200, n_ep = 2, arm.name = c("A","B"),
ep_type = "tte", cdist.rate = 0.5, sim_method = "copula",
copula_trt=copula::normalCopula(0.9), margins_trt=c("gamma", "beta"),
paramMargins_trt=list(list(shape=2, scale=1),list(shape1=2, shape2=2)),
copula_con=copula::normalCopula(0.9), margins_con=c("gamma", "beta"),
paramMargins_con=list(list(shape=2, scale=1),list(shape1=2, shape2=2)),
max_accrual_time = 5)

partition_t.plot(data, Ctime = seq(0,8,0.2), arm.name = c("A","B"),
priority = c(1,2), tau = 0, plotTimeUnit = "days", trt_group = "both")

#### An simulated example with three TTE endpoints.
data <- sim.data(n_trt = 200, n_con = 200, n_ep = 3, arm.name = c("A","B"),
ep_type = "tte", cdist.rate = 1, sim_method = "copula",
copula_trt=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3, dispstr = "un"),
margins_trt=c("gamma", "beta", "gamma"),
paramMargins_trt=list(list(shape=2, scale=2),list(shape1=2, shape2=2),list(shape=2, scale=3)),
copula_con=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3, dispstr = "un"),
margins_con=c("gamma", "beta", "gamma"),
paramMargins_con=list(list(shape=2, scale=1),list(shape1=2, shape2=1),list(shape=2, scale=2)),
max_accrual_time = 5)

partition_t.plot(data, Ctime = c(seq(0,8,0.5),seq(8.1,10,0.1)), arm.name = c("A","B"),
priority = c(3,2,1), tau = 0, plotTimeUnit = "years", trt_group = "trt")
}
}
\keyword{plot functions}
