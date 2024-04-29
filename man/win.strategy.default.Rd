\name{win.strategy.default}
\alias{win.strategy.default}
\title{
The Default Win Strategy Function.
}
\description{
An intermediate function to determine the win status for each pair based on the default win strategy. Specifically, one compares each subject in the treatment group with every subject in the control group to determine the win status.
}
\usage{
win.strategy.default(trt_con, priority, tau, np_direction)
}
\arguments{
  \item{trt_con}{Given N_t subjects in the treatment group and N_c subjects in the control group, there are N_t*N_c pairs, each row represents a pair. The analysis dataset trt_con contains the following variables:
    \describe{
      \item{stratum: }{A vector for the stratum number of the unmatched pairs.}
    \item{pid_trt: }{A vector for the subject id of the individuals from the treatment group within each unmatched pair.}
    \item{pid_con: }{A vector for the subject id of the individuals from the control group within each unmatched pair.}
    \item{Delta_j_trt: }{A vector for the event status of the j-th endpoint (1=event, 0=censored) for the individuals from the treatment group in each unmatched pair. If the outcome type for the endpoint is continuous/binary, then the event status is 1 for all.}
    \item{Delta_j_con: }{A vector for the event status of the j-th endpoint (1=event, 0=censored) for the individuals from the control group in each unmatched pair. If the outcome type for the endpoint is continuous/binary, then the event status is 1 for all.}
    \item{Y_j_trt: }{A vector for the outcome of the j-th endpoint for the individuals from the treatment group in each unmatched pair. For a time-to-event outcome, it would be a vector of observed time-to-event observations.}
    \item{Y_j_con: }{A vector for the outcome of the j-th endpoint for the individuals from the control group in each unmatched pair. For a time-to-event outcome, it would be a vector of observed time-to-event observations.}
    }
}
  \item{priority}{Importance order (from the most to the least important). For example, given three endpoints with the importance order as Endpoint 3, Endpoint 2, and Endpoint 1, input priority = c(3,2,1).}
  \item{tau}{A vector of numerical value for the magnitude of difference to determine win/loss/tie for each endpoint. Tau is applicable for TTE endpoints and continuous endpoints; tau is fixed as 0 for binary endpoints. Default is 0 for all endpoints.}
  \item{np_direction}{A vector of character for the direction to define a better result for each endpoint.
      \describe{
      \item{"larger": }{The treatment group wins if min(T_trt, C_trt, C_con + tau) > T_con + tau, for time-to-event outcomes; "the larger value wins" for continuous outcome, value 1 wins over the value 0 for binary outcome.}
      \item{"smaller": }{The treatment group wins if min(T_trt, C_trt, C_con - tau) < T_con - tau, for time-to-event outcomes; "the smaller value wins" for continuous outcome, value 0 wins over the value 1 for binary outcome.}
      }}
}
\details{
Users can define their own win strategy function and input the function with the argument "win.strategy" in the "win.stat" function.
}
\value{
\item{win_status}{A data frame for the win status of each pair for each endpoint.}
}
