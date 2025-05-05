\name{stat_t.plot}
\alias{stat_t.plot}
\title{
Plot Win Statistics over the Study Time.
}
\description{
Plot the win statistics as a function of the study time for time-to-event endpoints only. The function "win.stat" is called to calculate the win statistics at each time in the plot function.
}
\usage{
stat_t.plot(data, Ctime = Inf, arm.name = c(1,2), priority = c(1,2),
statistic = c("WR","NB","WO"),
Z_t_trt = NULL, Z_t_con = NULL, tau = 0,np_direction = "larger",horizon= Inf,
stratum.weight = c("unstratified","MH-type","wt.stratum1","wt.stratum2","equal"),
censoring_adjust = c("unadjusted","ipcw_tau","ipcw","covipcw"),
win.strategy = NULL, plotTimeUnit = NULL,
plot_CI = FALSE, alpha = 0.05, ...)
}
\arguments{
  \item{data}{The analysis dataset which contains the following variables:
    \describe{
      \item{arm: }{A vector for the treatment group (trt = 1 | 2), trt is the new treatment. Alternative names for "arm" include "trt", "treat" and "treatment".}
      \item{id: }{A vector for the patient id.}
      \item{stratum: }{A vector for the stratum.}
      \item{Delta_j: }{A vector for the event status of the j-th endpoint if the endpoint is a time-to-event outcome (1=event, 0=censored).}
      \item{Y_j: }{A vector for the outcome of the j-th endpoint, for a time-to-event outcome, it would be a vector for time.}
      \item{Start_time: }{A vector for the time when each of the individuals is first accrued to study.}
    }
  }
  \item{Ctime}{A vector of study times, at which the win proportions are to be plotted. Study time for a subject is the time since the subject is accrued, default as Inf.}
  \item{arm.name}{A vector for the labels of the two experimental arms, default to be c(1,2). The first label is for the treatment group, and the second label is for the control group.}
  \item{priority}{Importance order (from the most to the least important). For example, given three endpoints with the importance order as Endpoint 3, Endpoint 2, and Endpoint 1, input priority = c(3,2,1).}
  \item{statistic}{The win statistic to be plotted.
    \describe{
      \item{"WR": }{Win ratio.}
      \item{"NB": }{Net benefit.}
      \item{"WO": }{Win odds.}
    }
  }
  \item{Z_t_trt}{A matrix for the covariate history, each row is a (p+2) vector for one record of each subject in the treatment group. The first two columns are subject id, time, the other columns are the covariates (vector of length p). The baseline covariates are provided as the observed covariates corresponding to time 0.}
  \item{Z_t_con}{A matrix for the covariate history, each row is a (p+2) vector for one record of each subject in the control group. The first two columns are subject id, time, the other columns are the covariates (vector of length p). The baseline covariates are provided as the observed covariates corresponding to time 0.}
  \item{tau}{A vector of numerical value for the magnitude of difference to determine win/loss/tie for each endpoint. If tau is input as scalar, the function treat the taus for TTE endpoints and continuous endpoints to be the same and taus for binary endpoints as 0. Default as 0 for all endpoints.}
  \item{np_direction}{A vector of character for the direction to define a better result for each endpoint.
      \describe{
      \item{"larger": }{The treatment group wins if min(T_trt, C_trt, C_con + tau) > T_con + tau, for time-to-event outcomes; "the larger value wins" for continuous outcome, value 1 wins over the value 0 for binary outcome.}
      \item{"smaller": }{The treatment group wins if min(T_trt, C_trt, C_con - tau) < T_con - tau, for time-to-event outcomes; "the smaller value wins" for continuous outcome, value 0 wins over the value 1 for binary outcome.}
      }}
  \item{horizon}{Time horizon of interest.}
  \item{stratum.weight}{The weighting method for each stratum. Default is "unstratified" for unstratified analysis. A stratified analysis is performed if other weight option is specified. Other possible choices for this argument are listed below.
    \describe{
      \item{"MH-type": }{weight the wins with the reciprocal of the stratum size following the Mantel-Haenszel type stratified analysis as described in Dong et al. (2018).}
      \item{"wt.stratum1": }{weight the win statistics with weight equal to the number of subjects in each stratum divided by the total number of subjects.}
      \item{"wt.stratum2": }{weight the win statistics with weight equal to the number of subjects with events (of any TTE endpoint) in each stratum divided by the total number of subjects with events (of any TTE endpoint).}
      \item{"equal": }{set equal weights for all stratum.}
    }
  }
  \item{censoring_adjust}{The method to adjust censoring for the kernal functions. Possible choices are listed below.
    \describe{
      \item{"No": }{Without using the IPCW approach to dealing with the censoring.}
      \item{"IPCW": }{IPCW adjustment for censoring with the Kaplan-Meier estimator.}
      \item{"CovIPCW": }{IPCW adjustment for censoring with the time dependent Cox model.}
    }
  }
  \item{win.strategy}{The strategy to determine the win status. Default as NULL. If NULL, the default win strategy funtion "win.strategy.default" is called, see \link{win.strategy.default} for more details. Users can also define their own "win.strategy" function.}
  \item{plotTimeUnit}{The time units, e.g., "days", "months", "years", default as NULL.}
  \item{plot_CI}{If TRUE, plot the pointwise confidence interval, default as FALSE.}
  \item{alpha}{The significance level, default to be 0.05.}
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

stat_t.plot(data, arm.name = c("A","B"),priority = c(1,2),
Ctime = seq(2,12,1),plotTimeUnit = "years",statistic = "WR", tau = 0,
stratum.weight = "unstratified", plot_CI = TRUE)

#### An simulated example with three TTE endpoints.
data <- sim.data(n_trt = 200, n_con = 200, n_ep = 3, arm.name = c("A","B"),
ep_type = "tte", cdist.rate = 0.5, sim_method = "copula",
copula_trt=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3, dispstr = "un"),
margins_trt=c("gamma", "beta", "t"),
paramMargins_trt=list(list(shape=2, scale=1),list(shape1=2, shape2=2),list(df=5)),
copula_con=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3, dispstr = "un"),
margins_con=c("gamma", "beta", "t"),
paramMargins_con=list(list(shape=1, scale=1),list(shape1=2, shape2=3),list(df=5)),
max_accrual_time = 5)

stat_t.plot(data, arm.name = c("A","B"),priority = c(3,2,1),
Ctime = seq(1,8,0.5),plotTimeUnit = "years", statistic = "WR",
tau = 0, plot_CI = TRUE)
}
}
\keyword{plot functions}
