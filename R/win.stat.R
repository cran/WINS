win.stat<-function(data, ep_type, Z_t_trt = NULL, Z_t_con = NULL, arm.name = c(1,2),
                   priority = c(1,2), alpha = 0.05, digit = 5, tau = 0,
                   win.strategy = NULL,
                   pvalue = c("one-sided","two-sided"),
                   weight = c("unstratified","MH-type","wt.stratum1","wt.stratum2","equal"),
                   censoring_adjust = c("No","IPCW","CovIPCW"), var_method = c("Dong et al."),
                   summary.print = TRUE, ...){
  #### match the argument
  pvalue = match.arg(pvalue)
  weight = match.arg(weight)
  censoring_adjust = match.arg(censoring_adjust)

  #### obtain the number of endpoints and total number of individuals
  n_ep = length(priority)
  n_total = dim(data)[1]

  #### If ep_type is input as scalar, treat all the endpoints as the same type.
  if(length(ep_type)==1){
    cat("The outcome type for all the endpoints: ",ep_type,"\n")
    ep_type = rep(ep_type,n_ep)
  }else if(length(ep_type)!=n_ep){
    stop("The length of priority does not match the number of endpoints.")
  }

  #### If tau is input as scalar, treat all the tau as the same.
  if(length(tau)==1){
    tau = rep(tau,n_ep)
  }

  #############################################################################################
  #### Reorganize the data
  #############################################################################################
  colname.ds = colnames(data)

  if(max(c("arm","trt","treat","treatment")%in%colname.ds)==TRUE){
    arm = data[,which(colname.ds%in%c("arm","trt","treat","treatment"))]
  }else{
    stop("The treatment variable is not found. Please rename the treatment variable to arm, trt, treat or treatment.")
  }

  #### obtain the number of treatment and control
  n_trt = sum(arm==arm.name[1]); n_con = sum(arm==arm.name[2])

  if("id"%in%colname.ds==TRUE){
    id_trt = data[arm==arm.name[1],which(colname.ds%in%c("id"))]
    id_con = data[arm==arm.name[2],which(colname.ds%in%c("id"))]
  }else{
    id_trt = 1:n_trt; id_con = 1:n_con
  }

  if(("stratum"%in%colname.ds)==TRUE && weight != "unstratified"){
    stratum = data[,which(colname.ds=="stratum")]
  }else{
    # The unstratified win statistics are calculated as the special case for stratified analysis with only one stratum.
    stratum = rep(1,n_total)
    cat("This analysis is unstratified.","\n")
  }

  Delta = matrix(1,n_total,n_ep)
  if(max(c(stringr::str_detect(colname.ds,"Delta"),stringr::str_detect(colname.ds,"delta")))>0){
    ind_delta = which(stringr::str_detect(colname.ds,"Delta")|stringr::str_detect(colname.ds,"delta"))
    Delta[,which(ep_type=="tte")] = as.matrix(data[,ind_delta])
  }else if("tte"%in%ep_type){
    warning("No event status information detected. The default value of 1 is assigned to all individuals.")
  }

  Y = as.matrix(data[,which(stringr::str_detect(colname.ds,"Y"))])

  df = data.frame(arm = arm,stratum = stratum,Delta = Delta,Y = Y)
  colnames(df) = c("arm","stratum",paste0("Delta_",1:n_ep),paste0("Y_",1:n_ep))
  rm(arm,stratum,Delta,Y)

  #############################################################################################
  #### The options for data set with two endpoints and without stratum: obtain the win ratio
  #### and netbenefits from Luo et.al's method.
  #############################################################################################
  if(var_method == "Luo et al."){
    #### Obtain parameters for Luo et. al's method
    y1 = df[,which(colnames(df) == paste0("Y_",priority[2]))]
    y2 = df[,which(colnames(df) == paste0("Y_",priority[1]))]

    ind_delta_priority1 = which(colnames(df) == paste0("Delta_",priority[1]))
    ind_delta_priority2 = which(colnames(df) == paste0("Delta_",priority[2]))

    g = 1*(df$arm==arm.name[1])

    d = 1*(df[,ind_delta_priority1] == 1 & df[,ind_delta_priority2] == 1) +
        2*(df[,ind_delta_priority1] == 1 & df[,ind_delta_priority2] == 0) +
        3*(df[,ind_delta_priority1] == 0 & df[,ind_delta_priority2] == 0) +
        4*(df[,ind_delta_priority1] == 0 & df[,ind_delta_priority2] == 1)

    n = length(y1)

    #### Apply Luo et al.'s method
    luo_est = win.stat.Luo(y1,y2,d,g)

    wr = luo_est$wratio; sig2_log_wr = (n^(-1))*luo_est$vwratio/((luo_est$wratio)^2)
    nb = luo_est$wdiff/(n_trt*n_con); sig2_nb = (n^3)*luo_est$vwdiff/((n_trt*n_con)^2)

    #### p-value
    pvalue_WR = switch(pvalue,
                       "one-sided" = 1 - pnorm((log(wr)/sqrt(sig2_log_wr)), mean = 0, sd = 1),
                       "two-sided" = 2 - 2*pnorm(abs(log(wr)/sqrt(sig2_log_wr)), mean = 0, sd = 1))

    pvalue_NB = switch(pvalue,
                       "one-sided" = 1 - pnorm((nb/sqrt(sig2_nb)), mean = 0, sd = 1),
                       "two-sided" = 2 - 2*pnorm(abs(nb/sqrt(sig2_nb)), mean = 0, sd = 1))

    #### 100*(1-alpha)% CI
    WR_L = exp(log(wr) - qnorm(1-alpha/2)*sqrt(sig2_log_wr))
    WR_U = exp(log(wr) + qnorm(1-alpha/2)*sqrt(sig2_log_wr))

    NB_L = nb - qnorm(1-alpha/2)*sqrt(sig2_nb)
    NB_U = nb + qnorm(1-alpha/2)*sqrt(sig2_nb)

    Win_statisitc = list(Win_Ratio = list(WR = wr, WR_L = WR_L, WR_U = WR_U),
                         Net_Benefit = list(NB = nb, NB_L = NB_L, NB_U = NB_U))


    if(summary.print){
      #############################################################################################
      #### Output the result
      #############################################################################################
      cat("\n Win Ratio :", formatC(wr,digits = digit,format = "f"),"\n",
          "\n", pvalue, "p-value is: ",
          ifelse(pvalue_WR>10^(-digit),formatC(pvalue_WR,digits = digit,format = "f"),paste0("< ",10^(-digit))),
          "\n",
          "Lower limit of", 100*(1-alpha), "% CI of the win ratio: ", formatC(WR_L,digits = digit,format = "f"), "\n",
          "Upper limit of", 100*(1-alpha), "% CI of the win ratio: ", formatC(WR_U,digits = digit,format = "f"), "\n",
          "\n",
          "\n Net Benefit :", formatC(nb,digits = digit,format = "f"), "\n",
          "\n", pvalue, "p-value is: ",
          ifelse(pvalue_NB>10^(-digit),formatC(pvalue_NB,digits = digit,format = "f"),paste0("< ",10^(-digit))),
          "\n",
          "Lower limit of", 100*(1-alpha), "% CI of the net benefit: ", formatC(NB_L,digits = digit,format = "f"), "\n",
          "Upper limit of", 100*(1-alpha), "% CI of the net benefit: ", formatC(NB_U,digits = digit,format = "f"), "\n",
          "\n")
      return(invisible(NULL))
    }else{
      res = list(Win_statisitc = Win_statisitc, p_value = c(pvalue_WR,pvalue_NB))
      return(res)
    }
  }

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
    win_status = win.strategy(trt_con = trt_con, priority = priority, ...)
  }

  #############################################################################################
  #### Stratum-specific win ratios: estimated WR and the variance
  #############################################################################################
  #### Obtain kernel function K and L
  KL = switch (censoring_adjust,
               "No" = original.KL(win_status = win_status, trt_con = trt_con, n_ep = n_ep)$KL,
               "IPCW" = ipcw.adjusted.KL(win_status = win_status, trt = trt,con = con, trt_con = trt_con,
                                         priority = priority, n_ep = n_ep, ep_type = ep_type)$KL,
               "CovIPCW" = covipcw.adjusted.KL(win_status = win_status, trt = trt,con = con,trt_con = trt_con,
                                               Z_t_trt = Z_t_trt, Z_t_con = Z_t_con,
                                               priority = priority, n_ep = n_ep, ep_type = ep_type)$KL
  )

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
  w_stratum = switch(weight,
                     "unstratified" = 1,
                     "equal" = rep(1/length(N),length(N)),
                     "MH-type" = (1/N)/sum(1/N),
                     "wt.stratum1" = N/sum(N),
                     "wt.stratum2" = N_event/sum(N_event)
  )
  if(weight%in%c("unstratified","equal","MH-type")){
    stratified_N = sum((N_trt*N_con)*w_stratum)
  }
  stratified_WR = switch(weight,
                         "unstratified" = sum(win_trt*w_stratum/stratified_N)/sum(win_con*w_stratum/stratified_N),
                         "equal" = sum(win_trt*w_stratum/stratified_N)/sum(win_con*w_stratum/stratified_N),
                         "MH-type" = sum(win_trt*w_stratum/stratified_N)/sum(win_con*w_stratum/stratified_N),
                         "wt.stratum1" = sum(w_stratum*WR_stratum),
                         "wt.stratum2" = sum(w_stratum*WR_stratum)
  )
  stratified_NB = switch(weight,
                         "unstratified" = sum(win_trt*w_stratum/stratified_N)-sum(win_con*w_stratum/stratified_N),
                         "equal" = sum(win_trt*w_stratum/stratified_N)-sum(win_con*w_stratum/stratified_N),
                         "MH-type" = sum(win_trt*w_stratum/stratified_N)-sum(win_con*w_stratum/stratified_N),
                         "wt.stratum1" = sum(w_stratum*NB_stratum),
                         "wt.stratum2" = sum(w_stratum*NB_stratum)
  )
  stratified_WO = switch(weight,
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
  if(weight%in%c("wt.stratum1","wt.stratum2")){
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
    res = list(Win_prop = cbind(P_trt,P_con),Win_statistic = Win_statistic, p_value = c(pvalue_WR,pvalue_NB,pvalue_WO))
    return(res)
  }
}
