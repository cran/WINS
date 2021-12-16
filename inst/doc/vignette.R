## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=TRUE,echo=FALSE,message=FALSE,warning=FALSE-------------------------
library(WINS)

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
head(data_binary)

res_binary <- win.stat(data = data_binary, ep_type = "binary", priority = c(1:3),
                       weight = "unstratified", arm.name = c("A","B"), alpha=0.05, 
                       digit = 3,pvalue = "two-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
head(data_continuous)

res_continuous <- win.stat(data = data_continuous, ep_type = "continuous", 
                           arm.name = c("A","B"),tau = 0, priority = c(1:3), 
                           weight = "unstratified", alpha=0.05, digit = 3,
                           pvalue = "two-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
head(data_tte)

### Without using the IPCW approach to dealing with the censoring
res_tte <- win.stat(data = data_tte, ep_type = "tte", arm.name = c("A","B"),
                    tau = 0.1, priority = c(1:3), alpha = 0.05, digit = 3,
                    weight = "unstratified", censoring_adjust = "No", 
                    pvalue = "two-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
### IPCW adjustment for independent censoring
res_tte_ipcw <- win.stat(data = data_tte, ep_type = "tte", arm.name = c("A","B"),
                         tau = 0.1, priority = c(1:3), alpha = 0.05, digit = 3, 
                         weight = "unstratified", censoring_adjust = "IPCW", 
                         pvalue = "two-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
head(Z_t_trt)

### CovIPCW adjustment for dependent censoring
res_tte_covipcw <- win.stat(data = data_tte, ep_type = "tte", tau = 0.1, 
                            arm.name = c("A","B"), weight = "unstratified",
                            Z_t_trt = Z_t_trt, Z_t_con = Z_t_con,
                            priority = c(1:3), alpha = 0.05, digit = 3,
                            censoring_adjust = "CovIPCW", pvalue = "two-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
head(data_mix)

res_mix <- win.stat(data = data_mix, ep_type = c("tte","continuous","continuous"),
                    arm.name = c("A","B"), tau = 0.1, priority = c(1:3), 
                    alpha = 0.05, digit = 3, censoring_adjust = "No",
                    weight = "unstratified", pvalue = "two-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
### IPCW adjustment for independent censoring
res_mix_ipcw <- win.stat(data = data_mix, 
                         ep_type = c("tte","continuous","continuous"), 
                         arm.name = c("A","B"), tau = 0.1, priority = c(1:3),
                         alpha = 0.05, digit = 3, censoring_adjust = "IPCW",
                         weight="unstratified", pvalue = "two-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
res_mix_equal <- win.stat(data = data_mix_stratum, 
                          ep_type = c("tte","continuous","continuous"), 
                          arm.name = c("A","B"), tau = 0.1, priority = c(1:3), 
                          alpha = 0.05, digit = 3, censoring_adjust = "No", 
                          weight = "equal", pvalue = "one-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
head(data_mix_stratum)

res_mix_MH <- win.stat(data = data_mix_stratum, 
                       ep_type = c("tte","continuous","continuous"), 
                       arm.name = c("A","B"), tau = 0.1, priority = c(1:3), 
                       alpha = 0.05, digit = 3, censoring_adjust = "No", 
                       weight = "MH-type", pvalue = "one-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
res_mix_wt1 <- win.stat(data = data_mix_stratum, 
                        ep_type = c("tte","continuous","continuous"), 
                        arm.name = c("A","B"), tau = 0.1, priority = c(1:3), 
                        alpha = 0.05, digit = 3, censoring_adjust = "No", 
                        weight = "wt.stratum1", pvalue = "one-sided")

## ---- eval=TRUE,echo=TRUE,warning=FALSE---------------------------------------
res_mix_wt2 <- win.stat(data = data_mix_stratum, 
                        ep_type = c("tte","continuous","continuous"), 
                        arm.name = c("A","B"), tau = 0.1, priority = c(1:3), 
                        alpha = 0.05, digit = 3, censoring_adjust = "No", 
                        weight = "wt.stratum2", pvalue = "one-sided")

## ----eval=TRUE,echo=TRUE,warning=FALSE----------------------------------------
#### Generate two TTE endpoints with the more important TTE endpoint expected to occur
#### later with exponential distributions.
sim.data_tte <- sim.data(n_trt = 150, n_con = 100, n_ep = 2, 
                         arm.name = c("A","B"), ep_type = c("tte","tte"), 
                         cdist.rate = 0.5, sim_method = "tte_exponential", 
                         rate_trt = c(0.3,0.2), rate_con = c(0.5,0.4), 
                         max_accrual_time = 5)

win_stat_sim_tte <- win.stat(data = sim.data_tte, ep_type = c("tte","tte"), 
                             arm.name = c("A","B"), digit = 3, priority = c(2,1))

## ----eval=TRUE,echo=TRUE,warning=FALSE----------------------------------------
sim.data <- sim.data(n_trt = 150, n_con = 100, n_ep = 3, arm.name = c("A","B"),
                     ep_type = c("tte","tte","continuous"), 
                     cdist.rate = 0.5, sim_method = "copula",
                     copula_trt=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3,
                                             dispstr = "un"),
                     margins_trt=c("gamma", "beta", "t"),
                     paramMargins_trt=list(list(shape=2, scale=1),
                                           list(shape1=2, shape2=2), list(df=5)),
                     copula_con=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3,
                                             dispstr = "un"),
                     margins_con=c("gamma", "beta", "t"),
                     paramMargins_con=list(list(shape=1, scale=1),
                                           list(shape1=1, shape2=2), list(df=2)),
                     max_accrual_time = 5)

win_stat <- win.stat(data = sim.data, ep_type = c("tte","tte","continuous"),
                     arm.name = c("A","B"), digit = 3, priority = c(1,2,3))

## ----eval=TRUE,echo=TRUE,dpi=100,fig.height=4, fig.width=6,fig.align="center",message=FALSE,warning=FALSE----
#### An simulated example with three TTE endpoints.
data <- sim.data(n_trt = 200, n_con = 200, n_ep = 3, arm.name = c("A","B"),
                 ep_type = "tte", cdist.rate = 1, sim_method = "copula",
                 copula_trt=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3, 
                                         dispstr = "un"),
                 margins_trt=c("gamma", "beta", "gamma"),
                 paramMargins_trt=list(list(shape=2, scale=2),
                                       list(shape1=2, shape2=2),
                                       list(shape=2, scale=3)),
                 copula_con=copula::normalCopula(param=c(0.9,0.8,0.95), dim = 3, 
                                         dispstr = "un"),
                 margins_con=c("gamma", "beta", "gamma"),
                 paramMargins_con=list(list(shape=2, scale=1),
                                       list(shape1=2, shape2=1),
                                       list(shape=2, scale=2)),
                 max_accrual_time = 5)

stat_t.plot(data, arm.name = c("A","B"),priority = c(3,2,1), 
              Ctime = seq(2,8,0.2),plotTimeUnit = "years",
              statistic = "WR", tau = 0, plot_CI = TRUE)

## ----eval=TRUE,echo=TRUE,dpi=100,fig.height=4, fig.width=6,fig.align="center",message=FALSE,warning=FALSE----
partition_t.plot(data, Ctime = c(seq(0,9,0.5),seq(9.1,11,0.1)), arm.name = c("A","B"),
priority = c(3,2,1), tau = 0, plotTimeUnit = "years", trt_group = "trt")

