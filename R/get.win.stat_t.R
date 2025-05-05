get.win.stat_t<-function(trt, con, ep_type, Z_t_trt = NULL, Z_t_con = NULL, priority = c(1,2),
                         Ctimej = Inf, Start_time_trt = NULL, Start_time_con = NULL, alpha = 0.05,
                         tau = 0, np_direction = "larger",
                         stratum.weight = c("unstratified","MH-type","wt.stratum1","wt.stratum2","equal"),
                         censoring_adjust = c("unadjusted","ipcw_tau","ipcw","covipcw"),
                         pvalue = c("one-sided","two-sided"),
                         win.strategy = NULL, status_only = FALSE, return_CI = FALSE,
                         return_pvalue = FALSE, ...){
  #### obtain the number of endpoints and total number of individuals
  n_ep = length(priority)
  n_tte = sum(ep_type=="tte")

  #### Modify trt and con based on the follow time relative to the start of study
  #### and obtain the unmatched pairs
  colname.trt = colnames(trt)

  ind.delta1 = which(colname.trt == "Delta_1_trt")
  ind.time1 = which(colname.trt == "Y_1_trt")

  for(ind.ep in which(ep_type == "tte")){
    delta_new_trt = trt[,ind.delta1+ind.ep-1] * (trt[,ind.time1+ind.ep-1] <= (Ctimej - Start_time_trt))
    time_new_trt = apply(cbind(trt[,ind.time1+ind.ep-1],(Ctimej - Start_time_trt)),1,
                         func<-function(x) ifelse(min(x)>0,min(x),0))

    trt[,ind.delta1+ind.ep-1] = delta_new_trt
    trt[,ind.time1+ind.ep-1] = time_new_trt

    delta_new_con = con[,ind.delta1+ind.ep-1] * (con[,ind.time1+ind.ep-1] <= (Ctimej - Start_time_con))
    time_new_con = apply(cbind(con[,ind.time1+ind.ep-1],(Ctimej - Start_time_con)),1,
                         func<-function(x) ifelse(min(x)>0,min(x),0))

    con[,ind.delta1+ind.ep-1] = delta_new_con
    con[,ind.time1+ind.ep-1] = time_new_con
  }

  trt_con = merge(trt,con,by="stratum")

  #############################################################################################
  #### Determine winners/losers/ties
  #############################################################################################
  if(is.null(win.strategy)){
    win_status = win.strategy.default(trt_con = trt_con, priority = priority, tau = tau,
                                      np_direction = np_direction)
  }else{
    # user defined function
    win_status = win.strategy(trt_con = trt_con, priority = priority,
                              np_direction = np_direction, ...)
  }

  #### Obtain kernel function K and L
  res_KL = switch (censoring_adjust,
               "unadjusted" = original.KL(win_status = win_status, trt_con = trt_con, n_ep = n_ep),
               "ipcw_tau" = ipcw.adjusted.tau.KL(win_status = win_status, trt = trt, con = con, trt_con = trt_con,
                                         priority = priority, n_ep = n_ep,ep_type = ep_type,tau=tau),
               "ipcw" = ipcw.adjusted.KL(win_status = win_status, trt = trt, con = con, trt_con = trt_con,
                                         priority = priority, n_ep = n_ep, ep_type = ep_type),
               "covipcw" = covipcw.adjusted.KL(win_status = win_status, trt = trt,con = con,trt_con = trt_con,
                                               Z_t_trt = Z_t_trt, Z_t_con = Z_t_con,
                                               priority = priority, n_ep = n_ep, ep_type = ep_type)
  )

  KL = res_KL$KL
  status_KL = res_KL$status_KL

  #### number of patients per stratum
  N_trt = as.data.frame(table(trt$stratum))[,2]
  N_con = as.data.frame(table(con$stratum))[,2]

  N_trt_con = cbind(as.numeric(levels(factor(trt_con$stratum))), N_trt, N_con )
  colnames(N_trt_con)=c('stratum', 'N2_trt', 'N2_con')

  #### Estimate WR
  KL.summary = apply(cbind(KL$K,KL$L), 2, func<-function(x){
    temp = aggregate(x, by=list(Category=trt_con[[1]]), FUN=sum)
    as.matrix(temp)[,2]
  })
  KL.summary = matrix(KL.summary,ncol = 2)

  n_str = nrow(KL.summary)
  for(stri in 1:n_str){
    sum_KL_str = sum(KL.summary[stri,])/(N_trt[stri]*N_con[stri])
    if(sum_KL_str>=1){
      KL$K[which(KL$stratum==stri)] = KL$K[which(KL$stratum==stri)]/sum_KL_str
      KL$L[which(KL$stratum==stri)] = KL$L[which(KL$stratum==stri)]/sum_KL_str
      status_KL[which(KL$stratum==stri),] = status_KL[which(KL$stratum==stri),]/sum_KL_str
      KL.summary[stri,] = KL.summary[stri,]/sum_KL_str
    }
  }

  #### number of wins per stratum
  win_trt = KL.summary[,1]
  win_con = KL.summary[,2]

  if(status_only){
    return(list(win_status = status_KL))
  }else{
    #### stratum-specific win ratio, net benefit and win odds
    P_trt = win_trt/(N_trt*N_con)
    P_con = win_con/(N_trt*N_con)

    WR = win_trt/win_con
    NB = P_trt - P_con
    WO = (P_trt + 0.5*(1-P_trt-P_con))/(P_con + 0.5*(1-P_trt-P_con))

    #############################################################################################
    #### Stratified win ratio
    #############################################################################################

    #### Total sample size and total event per stratum
    N = N_trt + N_con

    ind.trt = which(colnames(trt)=="Delta_1_trt")
    event_trt = apply(as.matrix(trt[,c(ind.trt:(ind.trt+n_ep-1))],ncol = n_tte), 1, function(x) max(x)>0)
    N_event_trt = tapply(event_trt,trt$stratum,sum)

    ind.con = which(colnames(con)=="Delta_1_con")
    event_con = apply(as.matrix(con[,c(ind.con:(ind.con+n_ep-1))],ncol = n_tte), 1, function(x) max(x)>0)
    N_event_con = tapply(event_con,con$stratum,sum)

    N_event = N_event_trt + N_event_con

    #### Stratified win statistics
    w_stratum = switch(stratum.weight,
                       "unstratified" = 1,
                       "equal" = rep(1/length(N),length(N)),
                       "MH-type" = (1/N)/sum(1/N),
                       "wt.stratum1" = N/sum(N),
                       "wt.stratum2" = N_event/sum(N_event)
    )
    if(stratum.weight%in%c("unstratified","equal","MH-type")){
      stratified_N = sum((N_trt*N_con)*w_stratum)
    }
    stratified_WR = switch(stratum.weight,
                           "unstratified" = sum(win_trt*w_stratum/stratified_N)/sum(win_con*w_stratum/stratified_N),
                           "equal" = sum(win_trt*w_stratum/stratified_N)/sum(win_con*w_stratum/stratified_N),
                           "MH-type" = sum(win_trt*w_stratum/stratified_N)/sum(win_con*w_stratum/stratified_N),
                           "wt.stratum1" = sum(w_stratum*WR_stratum),
                           "wt.stratum2" = sum(w_stratum*WR_stratum)
    )
    stratified_NB = switch(stratum.weight,
                           "unstratified" = sum(win_trt*w_stratum/stratified_N)-sum(win_con*w_stratum/stratified_N),
                           "equal" = sum(win_trt*w_stratum/stratified_N)-sum(win_con*w_stratum/stratified_N),
                           "MH-type" = sum(win_trt*w_stratum/stratified_N)-sum(win_con*w_stratum/stratified_N),
                           "wt.stratum1" = sum(w_stratum*NB_stratum),
                           "wt.stratum2" = sum(w_stratum*NB_stratum)
    )
    stratified_WO = switch(stratum.weight,
                           "unstratified" = sum((sum(win_trt*w_stratum/stratified_N) +
                                                   0.5*(1-sum(win_trt*w_stratum/stratified_N)-
                                                          sum(win_con*w_stratum/stratified_N))))/
                             sum((sum(win_con*w_stratum/stratified_N) + 0.5*(1-sum(win_trt*w_stratum/stratified_N)-
                                                                               sum(win_con*w_stratum/stratified_N)))),
                           "equal" = sum((sum(win_trt*w_stratum/stratified_N) +
                                            0.5*(1-sum(win_trt*w_stratum/stratified_N)-
                                                   sum(win_con*w_stratum/stratified_N))))/
                             sum((sum(win_con*w_stratum/stratified_N) + 0.5*(1-sum(win_trt*w_stratum/stratified_N)-
                                                                               sum(win_con*w_stratum/stratified_N)))),
                           "MH-type" = sum((sum(win_trt*w_stratum/stratified_N) +
                                              0.5*(1-sum(win_trt*w_stratum/stratified_N)-
                                                     sum(win_con*w_stratum/stratified_N))))/
                             sum((sum(win_con*w_stratum/stratified_N) + 0.5*(1-sum(win_trt*w_stratum/stratified_N)-
                                                                               sum(win_con*w_stratum/stratified_N)))),
                           "wt.stratum1" = sum(w_stratum*WO_stratum),
                           "wt.stratum2" = sum(w_stratum*WO_stratum)
    )

    Win_statisitc = list(Win_Ratio = stratified_WR,
                         Net_Benefit = stratified_NB,
                         Win_Odds = stratified_WO,
                         win_status = win_status)

    if(return_CI == TRUE){
      #############################################################################################
      #### variances and covariances per stratum if there are any stratum.
      #############################################################################################
      #### calculate theta_K0/theta_L0
      theta_KL_0 = (win_trt + win_con)/(2*N_trt*N_con)
      theta_KL_0 = cbind(as.numeric(levels(factor(trt_con$stratum))), theta_KL_0)
      colnames(theta_KL_0)=c('stratum', 'theta_KL_0')

      sum_k_trt = aggregate(KL$K, by=list(trt_con$stratum, trt_con$pid_trt), FUN = sum)
      sum_k_con = aggregate(KL$K, by=list(trt_con$stratum, trt_con$pid_con), FUN = sum)

      sum_L_trt = aggregate(KL$L, by=list(trt_con$stratum, trt_con$pid_trt), FUN = sum)
      sum_L_con = aggregate(KL$L, by=list(trt_con$stratum, trt_con$pid_con), FUN = sum)

      names(sum_k_trt) = c('stratum', 'pid_trt', 'sum_k_trt')
      names(sum_k_con) = c('stratum', 'pid_con', 'sum_k_con')

      names(sum_L_trt) = c('stratum', 'pid_trt', 'sum_l_trt')
      names(sum_L_con) = c('stratum', 'pid_con', 'sum_l_con')

      KL = merge(KL, sum_k_trt, by=c('stratum', 'pid_trt'))
      KL = merge(KL, sum_k_con, by=c('stratum', 'pid_con'))
      KL = merge(KL, sum_L_trt, by=c('stratum', 'pid_trt'))
      KL = merge(KL, sum_L_con, by=c('stratum', 'pid_con'))
      KL = merge(KL, theta_KL_0, by=c('stratum'))
      KL = merge(KL, N_trt_con, by=c('stratum'))

      sig2_trt_1 = N_trt*N_con*aggregate( (KL$K-KL$theta_KL_0)*(KL$sum_k_trt - KL$K - (KL$N2_con - 1)*KL$theta_KL_0 ),
                                          by=list(KL$stratum), FUN = sum)[,2]/(N_con-1)
      sig2_trt_2 = N_trt*N_con*aggregate( (KL$K-KL$theta_KL_0)*(KL$sum_k_con - KL$K - (KL$N2_trt - 1)*KL$theta_KL_0 ),
                                          by=list(KL$stratum), FUN = sum)[,2]/(N_trt-1)

      sig2_con_1 = N_trt*N_con*aggregate( (KL$L-KL$theta_KL_0)*(KL$sum_l_con - KL$L - (KL$N2_trt - 1)*KL$theta_KL_0 ),
                                          by=list(KL$stratum), FUN = sum)[,2]/(N_trt-1)
      sig2_con_2 = N_trt*N_con*aggregate( (KL$L-KL$theta_KL_0)*(KL$sum_l_trt - KL$L - (KL$N2_con - 1)*KL$theta_KL_0 ),
                                          by=list(KL$stratum), FUN = sum)[,2]/(N_con-1)

      sig_trt_con_1 = N_trt*N_con*aggregate( (KL$K-KL$theta_KL_0)*(KL$sum_l_trt - KL$L - (KL$N2_con - 1)*KL$theta_KL_0 ),
                                             by=list(KL$stratum), FUN = sum)[,2]/(N_con-1)
      sig_trt_con_2 = N_trt*N_con*aggregate( (KL$K-KL$theta_KL_0)*(KL$sum_l_con - KL$L - (KL$N2_trt - 1)*KL$theta_KL_0 ),
                                             by=list(KL$stratum), FUN = sum)[,2]/(N_trt-1)

      sig2_trt = sig2_trt_1/N_trt + sig2_trt_2/N_con
      sig2_con = sig2_con_1/N_con + sig2_con_2/N_trt
      sig_trt_con = sig_trt_con_1/N_trt + sig_trt_con_2/N_con

      theta_tc = (win_trt + win_con)/2
      gam = theta_tc + 0.5*(N_trt*N_con-theta_tc-theta_tc)

      if(stratum.weight%in%c("wt.stratum1","wt.stratum2")){
        #### asymptotic variance
        sig2_wr = (sig2_trt + sig2_con - 2*sig_trt_con)/((theta_tc)^2)
        sig2_wo = (sig2_trt + sig2_con - 2*sig_trt_con)*(((N_trt*N_con)/(2*(gam^2)))^2)
        sig2_nb = (sig2_trt + sig2_con - 2*sig_trt_con)/((N_trt*N_con)^2)

        stratified_sig2_log_wr = sum(((w_stratum)^2)*sig2_wr)/(stratified_WR^2)
        stratified_sig2_log_wo = sum(((w_stratum)^2)*sig2_wo)/(stratified_WO^2)
        stratified_sig2_nb = sum(((w_stratum)^2)*sig2_nb)
      }else{
        #### asymptotic variance
        stratified_sig2_log_wr = sum(((w_stratum)^2)*(sig2_trt + sig2_con - 2*sig_trt_con))/
          ((sum(w_stratum*theta_tc))^2)
        stratified_sig2_log_wo = sum(((w_stratum)^2)*(sig2_trt + sig2_con - 2*sig_trt_con))*
          ((0.5/sum(w_stratum*gam)+0.5/sum(w_stratum*(N_trt*N_con-gam)))^2)
        stratified_sig2_nb = sum(((w_stratum)^2)*(sig2_trt + sig2_con - 2*sig_trt_con))/
          ((sum(w_stratum*N_trt*N_con))^2)
      }

      #### 100*(1-alpha)% CI
      stratified_WR_L = exp(log(stratified_WR) - qnorm(1-alpha/2)*sqrt(stratified_sig2_log_wr))
      stratified_WR_U = exp(log(stratified_WR) + qnorm(1-alpha/2)*sqrt(stratified_sig2_log_wr))

      stratified_WO_L = exp(log(stratified_WO) - qnorm(1-alpha/2)*sqrt(stratified_sig2_log_wo))
      stratified_WO_U = exp(log(stratified_WO) + qnorm(1-alpha/2)*sqrt(stratified_sig2_log_wo))

      stratified_NB_L = stratified_NB - qnorm(1-alpha/2)*sqrt(stratified_sig2_nb)
      stratified_NB_U = stratified_NB + qnorm(1-alpha/2)*sqrt(stratified_sig2_nb)

      CI_t = list(WR = c(stratified_WR_L,stratified_WR_U),
                  NB = c(stratified_NB_L,stratified_NB_U),
                  WO = c(stratified_WO_L,stratified_WO_U))

      if(return_pvalue){
        #### p-value
        pvalue_WR = switch(pvalue,
                           "one-sided" = 1 - pnorm((log(stratified_WR)/sqrt(stratified_sig2_log_wr)),
                                                   mean = 0, sd = 1),
                           "two-sided" = 2 - 2*pnorm(abs(log(stratified_WR)/sqrt(stratified_sig2_log_wr)),
                                                     mean = 0, sd = 1))

        pvalue_WO = switch(pvalue,
                           "one-sided" = 1 - pnorm((log(stratified_WO)/sqrt(stratified_sig2_log_wo)),
                                                   mean = 0, sd = 1),
                           "two-sided" = 2 - 2*pnorm(abs(log(stratified_WO)/sqrt(stratified_sig2_log_wo)),
                                                     mean = 0, sd = 1))

        pvalue_NB = switch(pvalue,
                           "one-sided" = 1 - pnorm((stratified_NB/sqrt(stratified_sig2_nb)),
                                                   mean = 0, sd = 1),
                           "two-sided" = 2 - 2*pnorm(abs(stratified_NB/sqrt(stratified_sig2_nb)),
                                                     mean = 0, sd = 1))
        p_value_t = c(pvalue_WR,pvalue_NB,pvalue_WO)
        return(list(Win_statisitc = Win_statisitc, p_value_t = p_value_t, CI_t = CI_t))
      }else{
        return(list(Win_statisitc = Win_statisitc, CI_t = CI_t))
      }
    }else{
      return(Win_statisitc = Win_statisitc)
    }
  }

}
