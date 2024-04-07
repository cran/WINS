unadjusted.win.stat<-function(df,n_total,arm.name = c(1,2),id_trt,id_con,
                              ep_type, n_ep,priority = c(1,2),tau = c(0,0),win.strategy = NULL,
                              alpha = 0.05,digit = 5,
                              pvalue, stratum.weight,
                              summary.print = TRUE, ...){
  #############################################################################################
  #### The options for general data set: obtain the estimate and variance estimate for
  #### win ratio, net benefit and win odds from Dong et.al's method.
  #############################################################################################
  #### pair the individuals in the treatment and control group

  trt = data.frame(id_trt,df[df$arm==arm.name[1],-1])
  colnames(trt) = c("pid_trt","stratum",paste0(colnames(df)[-c(1,2)],"_trt"))
  con = data.frame(id_con,df[df$arm==arm.name[2],-1])
  colnames(con) = c("pid_con","stratum",paste0(colnames(df)[-c(1,2)],"_con"))
  trt_con = merge(trt,con,by="stratum")

  #############################################################################################
  #### Determine winners/losers/ties
  #############################################################################################
  if(is.null(win.strategy)){
    win_status = win.strategy.default(trt_con = trt_con, priority = priority, tau = tau)
  }else{
    # user defined function
    win_status = win.strategy(trt_con = trt_con, priority = priority, tau=tau, ...)
  }

  #############################################################################################
  #### Stratum-specific win ratios: estimated WR and the variance
  #############################################################################################
  #### Obtain kernel function K and L
  KL = original.KL(win_status = win_status, trt_con = trt_con, n_ep = n_ep)$KL


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
      KL.summary[stri,] = KL.summary[stri,]/sum_KL_str
    }
  }

  #### number of wins per stratum
  win_trt = KL.summary[,1]
  win_con = KL.summary[,2]


  #### Obtain the number and the proportion of wins for each endpoint per stratum
  summary_ep = apply(win_status, 2, func<-function(x){
    temp1 = aggregate(x, by=list(trt_con$stratum), FUN = sum)
    temp2 = aggregate(x, by=list(trt_con$stratum), FUN = mean)
    temp_res = cbind(temp1,temp2[,2]); colnames(temp_res) = c("Stratum", "Count", "Proportion")
    return(temp_res)
  })


  #### stratum-specific win ratio, net benefit and win odds
  P_trt = win_trt/(N_trt*N_con)
  P_con = win_con/(N_trt*N_con)

  WR_stratum = win_trt/win_con
  NB_stratum = P_trt - P_con
  WO_stratum = (P_trt + 0.5*(1-P_trt-P_con))/(P_con + 0.5*(1-P_trt-P_con))

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
                                      by=list(KL$stratum), FUN = sum)[,2] / (N_con-1)
  sig2_trt_2 = N_trt*N_con*aggregate( (KL$K-KL$theta_KL_0)*(KL$sum_k_con - KL$K - (KL$N2_trt - 1)*KL$theta_KL_0 ),
                                      by=list(KL$stratum), FUN = sum)[,2] / (N_trt-1)

  sig2_con_1 = N_trt*N_con*aggregate( (KL$L-KL$theta_KL_0)*(KL$sum_l_con - KL$L - (KL$N2_trt - 1)*KL$theta_KL_0 ),
                                      by=list(KL$stratum), FUN = sum)[,2] / (N_trt-1)
  sig2_con_2 = N_trt*N_con*aggregate( (KL$L-KL$theta_KL_0)*(KL$sum_l_trt - KL$L - (KL$N2_con - 1)*KL$theta_KL_0 ),
                                      by=list(KL$stratum), FUN = sum)[,2] / (N_con-1)

  sig_trt_con_1 = N_trt*N_con*aggregate( (KL$K-KL$theta_KL_0)*(KL$sum_l_trt - KL$L - (KL$N2_con - 1)*KL$theta_KL_0 ),
                                         by=list(KL$stratum), FUN = sum)[,2] / (N_con-1)
  sig_trt_con_2 = N_trt*N_con*aggregate( (KL$K-KL$theta_KL_0)*(KL$sum_l_con - KL$L - (KL$N2_trt - 1)*KL$theta_KL_0 ),
                                         by=list(KL$stratum), FUN = sum)[,2] / (N_trt-1)

  sig2_trt = sig2_trt_1/N_trt + sig2_trt_2/N_con
  sig2_con = sig2_con_1/N_con + sig2_con_2/N_trt
  sig_trt_con = sig_trt_con_1/N_trt + sig_trt_con_2/N_con

  theta_tc = (win_trt + win_con)/2
  gam = theta_tc + 0.5*(N_trt*N_con-theta_tc-theta_tc)


  #############################################################################################
  #### Stratified win ratio
  #############################################################################################

  #### Total sample size and total event per stratum
  N = N_trt + N_con

  ind.tte = which(ep_type == "tte")
  n_tte = length(ind.tte)

  if(n_tte > 0){
    ind.trt = which(colnames(trt)=="Delta_1_trt")
    event_trt = apply(as.matrix(trt[,c(ind.trt:(ind.trt+n_ep-1))[ind.tte]],ncol = n_tte), 1, function(x) max(x)>0)
    N_event_trt = tapply(event_trt,trt$stratum,sum)

    ind.con = which(colnames(con)=="Delta_1_con")
    event_con = apply(as.matrix(con[,c(ind.con:(ind.con+n_ep-1))[ind.tte]],ncol = n_tte), 1, function(x) max(x)>0)
    N_event_con = tapply(event_con,con$stratum,sum)

    N_event = N_event_trt + N_event_con
  }else{
    N_event = N
  }

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

  #### Variance, CI and p-value
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

  #### z-statistic and p-value
  zstat_WR = log(stratified_WR)/sqrt(stratified_sig2_log_wr)
  pvalue_WR = switch(pvalue,
                     "one-sided" = 1 - pnorm(zstat_WR,
                                             mean = 0, sd = 1),
                     "two-sided" = 2 - 2*pnorm(abs(zstat_WR),
                                               mean = 0, sd = 1))

  zstat_WO = log(stratified_WO)/sqrt(stratified_sig2_log_wo)
  pvalue_WO = switch(pvalue,
                     "one-sided" = 1 - pnorm(zstat_WO,
                                             mean = 0, sd = 1),
                     "two-sided" = 2 - 2*pnorm(abs(zstat_WO),
                                               mean = 0, sd = 1))

  zstat_NB = stratified_NB/sqrt(stratified_sig2_nb)
  pvalue_NB = switch(pvalue,
                     "one-sided" = 1 - pnorm(zstat_NB,
                                             mean = 0, sd = 1),
                     "two-sided" = 2 - 2*pnorm(abs(zstat_NB),
                                               mean = 0, sd = 1))

  #### 100*(1-alpha)% CI
  stratified_WR_L = exp(log(stratified_WR) - qnorm(1-alpha/2)*sqrt(stratified_sig2_log_wr))
  stratified_WR_U = exp(log(stratified_WR) + qnorm(1-alpha/2)*sqrt(stratified_sig2_log_wr))

  stratified_WO_L = exp(log(stratified_WO) - qnorm(1-alpha/2)*sqrt(stratified_sig2_log_wo))
  stratified_WO_U = exp(log(stratified_WO) + qnorm(1-alpha/2)*sqrt(stratified_sig2_log_wo))

  stratified_NB_L = stratified_NB - qnorm(1-alpha/2)*sqrt(stratified_sig2_nb)
  stratified_NB_U = stratified_NB + qnorm(1-alpha/2)*sqrt(stratified_sig2_nb)

  Win_statistic = list(Win_Ratio = c(WR = stratified_WR,
                                     WR_L = stratified_WR_L,
                                     WR_U = stratified_WR_U),
                       Net_Benefit = c(NB = stratified_NB,
                                       NB_L = stratified_NB_L,
                                       NB_U = stratified_NB_U),
                       Win_Odds = c(WO = stratified_WO,
                                    WO_L = stratified_WO_L,
                                    WO_U = stratified_WO_U))

  if(summary.print){
    #############################################################################################
    #### Output the result
    #############################################################################################
    cat("\n Win Ratio :", formatC(stratified_WR,digits = digit,format = "f"),"\n",
        "\n", pvalue, "p-value is: ",
        ifelse(pvalue_WR>10^(-digit),formatC(pvalue_WR,digits = digit,format = "f"),paste0("< ",10^(-digit))),
        "\n",
        "Lower limit of", 100*(1-alpha), "% CI of the win ratio: ", formatC(stratified_WR_L,digits = digit, format = "f"), "\n",
        "Upper limit of", 100*(1-alpha), "% CI of the win ratio: ", formatC(stratified_WR_U,digits = digit, format = "f"), "\n",
        "\n",
        "\n Net Benefit :", formatC(stratified_NB,digits = digit,format = "f"), "\n",
        "\n", pvalue, "p-value is: ",
        ifelse(pvalue_NB>10^(-digit),formatC(pvalue_NB,digits = digit,format = "f"),paste0("< ",10^(-digit))),
        "\n",
        "Lower limit of", 100*(1-alpha), "% CI of the net benefit: ", formatC(stratified_NB_L,digits = digit, format = "f"), "\n",
        "Upper limit of", 100*(1-alpha), "% CI of the net benefit: ", formatC(stratified_NB_U,digits = digit, format = "f"), "\n",
        "\n",
        "\n Win Odds :", formatC(stratified_WO,digits = digit,format = "f"), "\n",
        "\n", pvalue, "p-value is: ",
        ifelse(pvalue_WO>10^(-digit),formatC(pvalue_WO,digits = digit,format = "f"),paste0("< ",10^(-digit))),
        "\n",
        "Lower limit of", 100*(1-alpha), "% CI of the win odds: ", formatC(stratified_WO_L,digits = digit, format = "f"), "\n",
        "Upper limit of", 100*(1-alpha), "% CI of the win odds: ", formatC(stratified_WO_U,digits = digit, format = "f"), "\n",
        "\n")
    return(invisible(NULL))
  }else{
    Win_prop = data.frame(stratum = as.numeric(levels(factor(trt_con$stratum))),P_trt=P_trt,P_con=P_con)
    z_statistic = data.frame(cbind(zstat_WR,zstat_NB,zstat_WO)); rownames(z_statistic) = "value"
    p_value = data.frame(cbind(pvalue_WR,pvalue_NB,pvalue_WO)); rownames(p_value) = "value"
    res = list(Win_prop = Win_prop,
               Win_statistic = Win_statistic,
               z_statistic = z_statistic,
               p_value = p_value,
               summary_ep = summary_ep)
    return(res)
  }
}
