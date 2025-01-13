\name{win.stat}
\alias{win.stat}
\title{
Function to Calculate the Win Statistics
}
\description{
Calculate the win statistics for a mixture type of outcomes including time-to-event outcome, continuous outcome and binary outcome.
}
\usage{
win.stat(data, ep_type, Z_t_trt = NULL, Z_t_con = NULL, iptw.weight = NULL,
arm.name = c(1,2),priority = c(1,2), alpha = 0.05, digit = 5, tau = 0,
np_direction = "larger", win.strategy = NULL,pvalue = c("one-sided","two-sided"),
stratum.weight = c("unstratified","MH-type","wt.stratum1","wt.stratum2","equal"),
method = c("unadjusted","ipcw","covipcw","iptw","iptw_ipcw","iptw_covipcw"),
summary.print = TRUE, ...)
}
\arguments{
  \item{data}{The analysis dataset which contains the following variables:
    \describe{
      \item{arm: }{A vector for the treatment groups. Alternative names for "arm" include "trt", "treat" and "treatment".}
      \item{id: }{A vector for the patient id.}
      \item{stratum: }{A vector for the stratum.}
      \item{Delta_j: }{A vector for the event status of the j-th endpoint if the endpoint is a time-to-event outcome (1=event, 0=censored).}
      \item{Y_j: }{A vector for the outcome of the j-th endpoint, for time-to-event outcome, Y_j is a vector for the observed time.}
      \item{Start_time: }{A vector for the time when each of the individuals is first accrued to study.}
    }
  }
  \item{ep_type}{A vector for the outcome type for each endpoint.  If scalar, the function will treat all the endpoints as the same type. The types of outcome include:
    \describe{
      \item{"tte": }{Time-to-event outcome, with the default win strategy: the treatment group wins if min(T_trt, C_trt, C_con + tau) > T_con + tau.}
      \item{"continuous": }{Continuous outcome, with the default win strategy: the treatment group wins if Y_trt > Y_con + tau.}
      \item{"binary": }{Binary outcome coded as 0/1, with the default win strategy: 1 is the winner over 0.}
    }
  }
  \item{Z_t_trt}{A matrix of the covariate history in the treatment group, each row is a (p+2) vector: the first two columns are subject id, time, the other p columns are the covariates (vector of length p). The baseline covariates are provided as the observed covariates corresponding to time 0.}
  \item{Z_t_con}{A matrix of the covariate history in the control group, each row is a (p+2) vector: the first two columns are subject id, time, the other p columns are the covariates (vector of length p). The baseline covariates are provided as the observed covariates corresponding to time 0.}
  \item{iptw.weight}{The weight assigned to each individual if method = "iptw".}
  \item{arm.name}{A vector for the labels of the two experimental arms, default to be c(1,2). The first label is for the treatment group, and the second label is for the control group.}
  \item{priority}{Importance order (from the most to the least important). For example, given three endpoints with the importance order as Endpoint 3, Endpoint 2, and Endpoint 1, input priority = c(3,2,1).}
  \item{alpha}{The significance level, default to be 0.05.}
  \item{digit}{The number of digits for the output, default to be 5.}
  \item{tau}{A vector of numerical value for the magnitude of difference to determine win/loss/tie for each endpoint. Tau is applicable for TTE endpoints and continuous endpoints; tau is fixed as 0 for binary endpoints. Default is 0 for all endpoints.}
  \item{np_direction}{A vector of character for the direction to define a better result for each endpoint.
      \describe{
      \item{"larger": }{The treatment group wins if min(T_trt, C_trt, C_con + tau) > T_con + tau, for time-to-event outcomes; "the larger value wins" for continuous outcome, value 1 wins over the value 0 for binary outcome.}
      \item{"smaller": }{The treatment group wins if min(T_trt, C_trt, C_con - tau) < T_con - tau, for time-to-event outcomes; "the smaller value wins" for continuous outcome, value 0 wins over the value 1 for binary outcome.}
      }}
  \item{win.strategy}{The strategy to determine the win status. Default as NULL. If NULL, the default win strategy funtion "win.strategy.default" is called, see \link{win.strategy.default} for more details. Users can also define their own "win.strategy" function.}
  \item{pvalue}{The p-value type: "one-sided" or "two-sided".}
  \item{stratum.weight}{The weighting method for each stratum. Default is "unstratified" for unstratified analysis. A stratified analysis is performed if other weight option is specified. Other possible choices for this argument are listed below.
    \describe{
      \item{"MH-type": }{weight the wins with the reciprocal of the stratum size following the Mantel-Haenszel type stratified analysis as described in Dong et al. (2018).}
      \item{"wt.stratum1": }{weight the win statistics with weight equal to the number of subjects in each stratum divided by the total number of subjects.}
      \item{"wt.stratum2": }{weight the win statistics with weight equal to the number of subjects with events (of any TTE endpoint) in each stratum divided by the total number of subjects with events (of any TTE endpoint).}
      \item{"equal": }{set equal weights for all stratum.}
    }
  }
  \item{method}{The method to adjust the kernal functions. Possible choices are listed below.
    \describe{
      \item{"unadjusted": }{Without using the IPCW approach to dealing with the censoring.}
      \item{"ipcw": }{IPCW adjustment for censoring with the Kaplan-Meier estimator.}
      \item{"covipcw": }{IPCW adjustment for censoring with the time dependent Cox model.}
      \item{"iptw": }{Inverse probability of treatment weighting (IPTW) method.}
      \item{"iptw_ipcw"}{Combine the IPCW adjustment for censoring with the Kaplan-Meier estimator
      and the Inverse probability of treatment weighting (IPTW) method.}
      \item{"iptw_covipcw"}{Combine the IPCW adjustment for censoring with the time dependent Cox
      model and the Inverse probability of treatment weighting (IPTW) method.}
    }
    Other methods may be added in future versions.
  }
  \item{summary.print}{If TRUE, print out a summary of the estimation and inference result for the win statistics; If FALSE, return a list that summarizes the results. Default as TRUE.}
  \item{...}{Argument passed from user defined functions "win.strategy" if there is any. For instructions on this "win.strategy" function, see \link{win.strategy.default} for more details.}
}
\details{
The arguments of user defined "win.strategy" function must at least include the argument "trt_con" and "priority". "priority" is defined the same as stated in the main function "win.stat". The intermediate analysis dataset "trt_con" for the patient pairs (i.e., unmatched pairs, see Pocock et al., 2012) contains the following variables. Each row represents a pair.
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
\value{
  \item{Win_prop}{The win proportion of the treatment and the control group.}
  \item{Win_statistic}{The win statistics including:
    \describe{
      \item{Win_Ratio: }{A list for the ratio of the win proportion between the treatment and the control groups.}
      \item{Net_Benefit: }{A list for the difference of the win proportion between the treatment and the control groups.}
      \item{Win_Odds: }{A list for the win odds between the treatment and the control groups.}
    }
  }
  \item{z_statistic}{The z-scores including:
    \describe{
      \item{zstat_WR: }{z-score for win ratio.}
      \item{zstat_NB: }{z-score for net benefit.}
      \item{zstat_WO: }{z-score for win odds.}
    }
  }
  \item{pvalue}{The p-value for the test statistics including:
    \describe{
      \item{pvalue_WR: }{p-value for win ratio.}
      \item{pvalue_NB: }{p-value for net benefit.}
      \item{pvalue_WO: }{p-value for win odds.}
    }
  }
  \item{summary_ep}{The win count and win proportion of the treatment and the control group for each endpoint}
}
\references{
  \itemize{
    \item Bebu, I. and Lachin, J.M., 2016. Large sample inference for a win ratio analysis of a composite outcome based on prioritized components. Biostatistics.
    \item Buyse, M., 2010. Generalized pairwise comparisons of prioritized outcomes in the two-sample problem. Statistics in medicine.
    \item Cui, Y., Dong, G., Kuan, P.F. and Huang, B., 2022. Evidence synthesis analysis with prioritized benefit outcomes in oncology clinical trials. Journal of Biopharmaceutical Statistics.
    \item Dong, G., Li, D., Ballerstedt, S. and Vandemeulebroecke, M., 2016. A generalized analytic solution to the win ratio to analyze a composite endpoint considering the clinical importance order among components. Pharmaceutical statistics.
    \item Dong, G., Qiu, J., Wang, D. and Vandemeulebroecke, M., 2018. The stratified win ratio. Journal of biopharmaceutical statistics.
    \item Dong, G., Huang, B., Chang, Y.W., Seifu, Y., Song, J. and Hoaglin, D.C., 2020. The win ratio: Impact of censoring and follow-up time and use with nonproportional hazards. Pharmaceutical statistics.
    \item Dong, G., Mao, L., Huang, B., Gamalo-Siebers, M., Wang, J., Yu, G. and Hoaglin, D.C., 2020. The inverse-probability-of-censoring weighting (IPCW) adjusted win ratio statistic: an unbiased estimator in the presence of independent censoring. Journal of biopharmaceutical statistics.
    \item Dong, G., Hoaglin, DC., Qiu, J., Matsouaka, RA., Chang, Y., Wang, J., Vandemeulebroecke, M., 2020. The win ratio: on interpretation and handling of ties. Statistics in Biopharmaceutical Research.
    \item Dong, G., Huang, B., Wang, D., Verbeeck, J., Wang, J. and Hoaglin, D.C., 2021. Adjusting win statistics for dependent censoring. Pharmaceutical Statistics.
    \item Dong, G., Huang, B., Verbeeck, J., Cui, Y., Song, J., Gamalo-Siebers, M., Wang, D., Hoaglin, D.C., Seifu, Y., Mutze, T. and Kolassa, J., 2022. Win statistics (win ratio, win odds, and net benefit) can complement one another to show the strength of the treatment effect on time-to-event outcomes. Pharmaceutical Statistics.
    \item Dong, G., Hoaglin, D.C., Huang, B., Cui, Y., Wang, D., Cheng, Y. and Gamalo-Siebers, M., 2023. The stratified win statistics (win ratio, win odds, and net benefit). Pharmaceutical Statistics.
    \item Finkelstein, D.M. and Schoenfeld, D.A., 1999. Combining mortality and longitudinal measures in clinical trials. Statistics in medicine.
    \item Finkelstein, D.M. and Schoenfeld, D.A., 2019. Graphing the Win Ratio and its components over time. Statistics in medicine.
    \item Luo, X., Tian, H., Mohanty, S. and Tsai, W.Y., 2015. An alternative approach to confidence interval estimation for the win ratio statistic. Biometrics.
    \item Pocock, S.J., Ariti, C.A., Collier, T.J. and Wang, D., 2012. The win ratio: a new approach to the analysis of composite endpoints in clinical trials based on clinical priorities. European heart journal.
    \item Wang D, Pocock S., 2016. A win ratio approach to comparing continuous non-normal outcomes in clinical trials. Pharmaceutical Statistics.
    \item Wang, D., Zheng S., Cui, Y., He, N., Chen, T., Huang, B., 2023. Adjusted win ratio using inverse probability treatment weighting (IPTW) propensity score analysis. Journal of Biopharmaceutical Statistics.
    \item Zheng S, Wang D, Qiu J, Chen T, Gamalo M., 2023. A win ratio approach for comparing crossing survival curves in clinical trials. Journal of Biopharmaceutical Statistics.
  }
}


\examples{
\donttest{
#### An example with time-to-event outcome.
head(data_tte)
head(Z_t_trt)

### CovIPCW adjustment for dependent censoring
res_tte_covipcw <- win.stat(data = data_tte, ep_type = "tte", arm.name = c("A","B"), tau = 0.1,
Z_t_trt = Z_t_trt, Z_t_con = Z_t_con, priority = c(1:3), alpha = 0.05, digit = 3,
method = "covipcw", stratum.weight = "unstratified", pvalue = "two-sided")

#### An example with continuous outcome.
head(data_continuous)

res_continuous <- win.stat(data = data_continuous, ep_type = "continuous", arm.name = c("A","B"),
tau = 0, priority = c(1:3), alpha=0.05, digit = 3, stratum.weight = "unstratified",
pvalue = "two-sided")

#### An example with binary outcome.
head(data_binary)

res_binary <- win.stat(data = data_binary, ep_type = "binary", arm.name = c("A","B"),
priority = c(1:3), alpha=0.05, digit = 3, stratum.weight = "unstratified",
pvalue = "two-sided")
}
}
