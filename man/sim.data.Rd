\name{sim.data}
\alias{sim.data}
\title{
Function for Data Simulation
}
\description{
A function to simulate a dataset with dependent endpoints. The time-to-event endpoints generated are assumed to have noninformative censoring.
}
\usage{
sim.data(randomseed = 12345, n_trt = 200, n_con = 200, n_ep = 2, n_stratum = 1,
arm.name = c(1,2), ep_type, cdist.rate, sim_method = "copula",
copula_trt = NULL, margins_trt = NULL, paramMargins_trt = NULL,
copula_con = NULL, margins_con = NULL, paramMargins_con = NULL,
rate_trt = NULL, rate_con = NULL, max_accrual_time = NULL)
}
\arguments{
  \item{randomseed}{The random seed.}
  \item{n_trt}{The number of individuals in the treatment group.}
  \item{n_con}{The number of individuals in the control group.}
  \item{n_ep}{The number of endpoints.}
  \item{n_stratum}{The number of strata. For the simulated dataset, n_stratum is fixed at 1 assuming homogeneous population.}
  \item{arm.name}{A vector for the labels of the two experimental arms, default to be c(1,2). The first label is for the treatment group, and the second label is for the control group.}
  \item{ep_type}{A vector for the outcome type for each endpoint.  If scalar, the function will treat all the endpoints as the same type. The types of outcome include:
    \describe{
      \item{"tte": }{Time-to-event outcome, with the default win strategy: the treatment group wins if min(T_trt, C_trt, C_con + tau) > T_con + tau.}
      \item{"continuous": }{Continuous outcome, with the default win strategy: the treatment group wins if Y_trt > Y_con + tau.}
      \item{"binary": }{Binary outcome coded as 0/1, with the default win strategy: 1 is the winner over 0.}
    }
  }
  \item{cdist.rate}{The censoring time is generated from an exponential distribution. This argument is a vector with the rate of the censoring distribution for each time-to-event endpoint. If scalar, the function will treat all the rate for censoring distribution as the same.}
  \item{sim_method}{Method used to generate multivariate dependence. Possible choices include "copula" and "tte_exponential"}
  \item{copula_trt}{an object of "copula" for the treatment group.}
  \item{margins_trt}{a character vector specifying all the parametric marginal distributions for the treatment group. See details in the R documentation for function "copula::Mvd".}
  \item{paramMargins_trt}{a list for which each element is a list (or numeric vectors) of named components, giving the parameter values of the marginal distributions for the treatment group. See details in the R documentation for function "copula::Mvd".}
  \item{copula_con}{Same argument as "copula_trt" for the control group.}
  \item{margins_con}{Same argument as "margins_trt" for the control group.}
  \item{paramMargins_con}{Same argument as "paramMargins_trt" for the control group.}
  \item{rate_trt}{A vector of the rate in the treatment group for each time-to-event endpoint following an exponential distribution when "sim_method" is set to be the option "tte_exponential".}
  \item{rate_con}{A vector of the rate in the control group for each time-to-event endpoint following an exponential distribution when "sim_method" is set to be the option "tte_exponential".}
  \item{max_accrual_time}{if specified, simulate the study entry time for each individual from uniform distribution U(0,max_accrual_time).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
To learn more about "copula", please refer to a discussion on modelling dependence with copulas with the link \url{https://datascienceplus.com/modelling-dependence-with-copulas/}. It shows on a high level how copula works, how to use a copula in R using the copula package and then provides a simple example. Moreover, when "sim_method" is set to be the option "tte_exponential", we simulate two endpoints based on the exponential distribution. Dependence between the two simulated endpoints is introduced, as the earlier endpoint takes the min of the two simulated exponential variables.
}
\value{
  \item{data}{The analysis dataset which contains the following variables:
    \describe{
      \item{arm: }{A vector for the treatment group (trt = 1 | 2), 1 represents the treatment group and 2 represents the control group.}
      \item{stratum: }{A vector for the stratum number. Alternative names for "stratum" include "group", "level" and "grade".}
      \item{Delta_j: }{A vector for the event status of the j-th endpoint if the endpoint is time-to-event outcome (1=event, 0=censored).}
      \item{Y_j: }{A vector for the outcome of the j-th endpoint, for time-to-event outcome, it would be a vector of simulated time.}
      \item{Start_time: }{A vector for the time when each of the individuals is first accrued to study. Valid only if "max_accrual_time" is not NULL.}
    }
  }
}
\examples{
\donttest{
#### Generate with copula: This example is for three endpoints, noted as Y_1, Y_2, and Y_3,
#### with endpoint type as TTE, TTE and continuous.
#### For both the treatment group and the control group, the correlation coefficients
#### cor(Y_1,Y_2), cor(Y_1,Y_3) and cor(Y_2,Y_3) are 0.9, 0.8 and 0.95, respectively.
#### For each treatment group, the marginal distribution for Y_1, Y_2, and Y_3 are Gamma,
#### Beta and Student t specified as a vector in "margins_trt"/"margins_con". The parameters
#### are specified as a list corresponding to the margianl distributions in "paramMargins_trt"
#### or "paramMargins_con".
sim.data <- sim.data(n_trt = 150, n_con = 100, n_ep = 3, arm.name = c("A","B"),
ep_type = c("tte","tte","continuous"), cdist.rate = 0.5, sim_method = "copula",
copula_trt=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3, dispstr = "un"),
margins_trt=c("gamma", "beta", "t"),
paramMargins_trt=list(list(shape=2, scale=1),list(shape1=2, shape2=2),list(df=5)),
copula_con=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3, dispstr = "un"),
margins_con=c("gamma", "beta", "t"),
paramMargins_con=list(list(shape=1, scale=1),list(shape1=1, shape2=2),list(df=2)),
max_accrual_time = 5)

win_stat <- win.stat(data = sim.data, ep_type = c("tte","tte","continuous"),
arm.name = c("A","B"), priority = c(1,2,3))

#### Generate two TTE endpoints with the more important TTE endpoint expected to occur later
#### with exponential distribution.
sim.data2 <- sim.data(n_trt = 150, n_con = 100, n_ep = 2, arm.name = c("A","B"),
ep_type = c("tte","tte"), cdist.rate = 0.5, sim_method = "tte_exponential",
rate_trt = c(0.2,0.25),rate_con = c(0.4,0.5), max_accrual_time = 5)

win_stat2 <- win.stat(data = sim.data2, ep_type = c("tte","tte"), arm.name = c("A","B"),
priority = c(1,2))
}
}
