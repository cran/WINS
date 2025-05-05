ipcw.adjusted.tau.KL<-function(win_status, trt, con, trt_con, priority, n_ep, ep_type, tau){
  #### Obtain the number of treatment and control
  num1 = dim(trt)[1]
  num2 = dim(con)[1]

  #### Calculate K and L
  status_KL = matrix(0,nrow = nrow(win_status),ncol = ncol(win_status))

  # Setup survival object
  ind.tte = which(ep_type=="tte"); ntte = length(ind.tte)
  if(ntte != length(ep_type)){
    stop("The IPCW-adjusted approach can only work for the case where all endpoints are tte!")
  }
  cen_time_trt = apply(as.matrix(trt[,(2+n_ep+ind.tte)],ncol = ntte),1,max)
  cen_time_con = apply(as.matrix(con[,(2+n_ep+ind.tte)],ncol = ntte),1,max)

  event_trt = apply(as.matrix(trt[,(2+ind.tte)],ncol = ntte),1,min)
  event_con = apply(as.matrix(con[,(2+ind.tte)],ncol = ntte),1,min)

  trt_id = as.factor(trt_con$pid_trt); levels(trt_id) = 1:num1
  con_id = as.factor(trt_con$pid_con); levels(con_id) = 1:num2
  cen_time_trt_long = cen_time_trt[trt_id]
  cen_time_con_long = cen_time_con[con_id]
  event_trt_long = event_trt[trt_id]
  event_con_long = event_con[con_id]

  ind_trt = seq(1,num1*num2,num1)
  ind_con = 1:num2

  # obtain the compare matrix and vector
  compare.mat_trt = (cen_time_trt>=t(array(rep(cen_time_trt, num1), c(num1, num1))))*1
  compare.vec_trt = apply(compare.mat_trt,2,sum)

  compare.mat_con = (cen_time_con>=t(array(rep(cen_time_con, num2), c(num2, num2))))*1
  compare.vec_con = apply(compare.mat_con,2,sum)

  # Obtain Kaplan-Meier estimator, the survival rate is truncated by 0.1 if small in tail.
  surv_trt = survival::Surv(time = cen_time_trt, event = as.numeric(event_trt==0))
  csurv1.fit = survival::survfit(surv_trt~1,type='kaplan-meier')
  ind.trun1 = (csurv1.fit$surv >= 0.1)
  csurv1.fit$time = c(0,csurv1.fit$time[ind.trun1])
  csurv1.fit$surv = c(1,csurv1.fit$surv[ind.trun1])
  L.tpt_trt = length(csurv1.fit$time)

  surv_con = survival::Surv(time = cen_time_con, event = as.numeric(event_con==0))
  csurv2.fit = survival::survfit(surv_con~1,type='kaplan-meier')
  ind.trun2 = (csurv2.fit$surv >= 0.1)
  csurv2.fit$time = c(0,csurv2.fit$time[ind.trun2])
  csurv2.fit$surv = c(1,csurv2.fit$surv[ind.trun2])
  L.tpt_con = length(csurv2.fit$time)

  tte.count = 0; weight_t = weight_c = matrix(rep(1,nrow(trt_con)),ncol=1); signal = 1; pm_ind = NULL
  csurv1.x_trt = csurv1.x_con = csurv1.x_con_ptau = csurv1.x_con_mtau = NULL
  csurv2.x_trt = csurv2.x_trt_ptau = csurv2.x_trt_mtau = csurv2.x_con = NULL
  cen_time_trt_track = NULL; cen_time_con_track = NULL
  for(l in 1:n_ep){
    ind.prior = priority[l]; tau_l = tau[l]
    if(l != 1){
      lcount = 2^(l-2)
      weight_t <- cbind(weight_t,weight_t)
      weight_t[!((cen_time_trt_last>=cen_time_con_last-tau_l)*event_con_last == 1),1:lcount] = 0
      weight_t[!((cen_time_trt_last>cen_time_con_last+tau_l)*event_con_last == 1),(lcount+1):(2*lcount)] = 0


      weight_c <- cbind(weight_c,weight_c)
      weight_c[!((cen_time_con_last>=cen_time_trt_last-tau_l)*event_trt_last == 1),1:lcount] = 0
      weight_c[!((cen_time_con_last>cen_time_trt_last+tau_l)*event_trt_last == 1),(lcount+1):(2*lcount)] = 0


      signal <- c(signal,-signal)
      pm_ind <- paste0(rep(c("m","p"),each=lcount),rep(pm_ind,2))
    }
    # Setup survival object
    cen_time_trt_epl = trt[,(2+n_ep+ind.prior)]
    cen_time_con_epl = con[,(2+n_ep+ind.prior)]

    cen_time_trt_last <- trt_con[,(2+n_ep+ind.prior)]
    cen_time_con_last <- trt_con[,(3+3*n_ep+ind.prior)]

    event_trt_last <- trt_con[,(2+ind.prior)]
    event_con_last <- trt_con[,(3+2*n_ep+ind.prior)]

    cen_time_trt_track <- cbind(cen_time_trt_track,trt_con[,(2+n_ep+ind.prior)])
    cen_time_con_track <- cbind(cen_time_con_track,trt_con[,(3+3*n_ep+ind.prior)])

    csurv1.indx_trt = apply(cen_time_trt_epl >= t(array(rep(csurv1.fit$time, num1), c(L.tpt_trt, num1))), 1, sum)
    csurv1.x_trt = matrix(cbind(csurv1.x_trt,csurv1.fit$surv[csurv1.indx_trt]),ncol = l)

    csurv1.indx_con = apply(cen_time_con_epl >= t(array(rep(csurv1.fit$time, num2), c(L.tpt_trt, num2))), 1, sum)
    csurv1.x_con = matrix(cbind(csurv1.x_con,csurv1.fit$surv[csurv1.indx_con]),ncol = l)

    csurv1.indx_con_ptau = apply(cen_time_con_epl + tau_l >= t(array(rep(csurv1.fit$time, num2), c(L.tpt_trt, num2))), 1, sum)
    csurv1.indx_con_ptau = ifelse(csurv1.indx_con_ptau == 0, 1, csurv1.indx_con_ptau)
    csurv1.x_con_ptau = matrix(cbind(csurv1.x_con_ptau,
                                      csurv1.fit$surv[csurv1.indx_con_ptau]),ncol = l)

    csurv1.indx_con_mtau = apply(cen_time_con_epl - tau_l >= t(array(rep(csurv1.fit$time, num2), c(L.tpt_trt, num2))), 1, sum)
    csurv1.indx_con_mtau = ifelse(csurv1.indx_con_mtau == 0, 1, csurv1.indx_con_mtau)
    csurv1.x_con_mtau = matrix(cbind(csurv1.x_con_mtau,
                                      csurv1.fit$surv[csurv1.indx_con_mtau]),ncol = l)

    csurv2.indx_trt = apply(cen_time_trt_epl >= t(array(rep(csurv2.fit$time, num1), c(L.tpt_con, num1))), 1, sum)
    csurv2.x_trt = matrix(cbind(csurv2.x_trt,csurv2.fit$surv[csurv2.indx_trt]),ncol = l)

    csurv2.indx_con = apply(cen_time_con_epl >= t(array(rep(csurv2.fit$time, num2), c(L.tpt_con, num2))), 1, sum)
    csurv2.x_con = matrix(cbind(csurv2.x_con,csurv2.fit$surv[csurv2.indx_con]),ncol = l)

    csurv2.indx_trt_ptau = apply(cen_time_trt_epl + tau_l >= t(array(rep(csurv2.fit$time, num1), c(L.tpt_con, num1))), 1, sum)
    csurv2.indx_trt_ptau = ifelse(csurv2.indx_trt_ptau == 0, 1, csurv2.indx_trt_ptau)
    csurv2.x_trt_ptau = matrix(cbind(csurv2.x_trt_ptau,
                                      csurv2.fit$surv[csurv2.indx_trt_ptau]),ncol = l)

    csurv2.indx_trt_mtau = apply(cen_time_trt_epl - tau_l >= t(array(rep(csurv2.fit$time, num1), c(L.tpt_con, num1))), 1, sum)
    csurv2.indx_trt_mtau = ifelse(csurv2.indx_trt_mtau == 0, 1, csurv2.indx_trt_mtau)
    csurv2.x_trt_mtau = matrix(cbind(csurv2.x_trt_mtau,
                                      csurv2.fit$surv[csurv2.indx_trt_mtau]),ncol = l)

    if(l == 1){
      IPCW_t = csurv1.x_con_ptau[con_id]*csurv2.x_con[con_id]
      IPCW_c = csurv1.x_trt[trt_id]*csurv2.x_trt_ptau[trt_id]

      I_K_ijv = weight_t*win_status[,(2*l-1)]
      g1Xc_K = cen_time_con_track + tau_l
      g2Xc_K = cen_time_con_track

      I_L_ijv = weight_c*win_status[,(2*l)]
      g1Xt_L = cen_time_trt_track
      g2Xt_L = cen_time_trt_track + tau_l

      status_KL[,(2*l-1)] = get.Kijv(num1,num2,trt_id, con_id, ind_trt, ind_con, IPCW_t, I_K_ijv, g1Xc_K, g2Xc_K,
                                     cen_time_trt, cen_time_con, cen_time_trt_long, cen_time_con_long,
                                     event_trt, event_con, event_trt_long, event_con_long,
                                     compare.vec_trt, compare.vec_con)
      status_KL[,(2*l)] = get.Lijv(num1,num2,trt_id, con_id, ind_trt, ind_con, IPCW_c, I_L_ijv, g1Xt_L, g2Xt_L,
                                   cen_time_trt, cen_time_con, cen_time_trt_long, cen_time_con_long,
                                   event_trt, event_con, event_trt_long, event_con_long,
                                   compare.vec_trt, compare.vec_con)
    }else{
      pm_ind_split = strsplit(pm_ind,split = "")
      for(part in 1:(2^(l-1))){
        part_pm = pm_ind_split[[part]]
        IPCW_t_matrix1 = csurv1.x_con_ptau[con_id,l]; IPCW_c_matrix2 = csurv2.x_trt_ptau[trt_id,l]
        IPCW_t_matrix2 = csurv2.x_con[con_id,]; IPCW_c_matrix1 = csurv1.x_trt[trt_id,]
        Xt_matrix = cen_time_trt_track[,l] + tau_l
        Xc_matrix = cen_time_con_track[,l] + tau_l
        for(pmi in 1:(l-1)){
          if(part_pm[pmi] == "m"){
            IPCW_t_matrix1 = cbind(IPCW_t_matrix1,csurv1.x_con_mtau[con_id,l-pmi])
            IPCW_c_matrix2 = cbind(IPCW_c_matrix2,csurv2.x_trt_mtau[trt_id,l-pmi])
            Xt_matrix = cbind(Xt_matrix,cen_time_trt_track[,l-pmi] - tau_l)
            Xc_matrix = cbind(Xc_matrix,cen_time_con_track[,l-pmi] - tau_l)
          }else{
            IPCW_t_matrix1 = cbind(IPCW_t_matrix1,csurv1.x_con_ptau[con_id,l-pmi])
            IPCW_c_matrix2 = cbind(IPCW_c_matrix2,csurv2.x_trt_ptau[trt_id,l-pmi])
            Xt_matrix = cbind(Xt_matrix,cen_time_trt_track[,l-pmi] + tau_l)
            Xc_matrix = cbind(Xc_matrix,cen_time_con_track[,l-pmi] + tau_l)
          }
        }
        IPCW_t = apply(IPCW_t_matrix1,1,min)*apply(IPCW_t_matrix2,1,min)
        IPCW_c = apply(IPCW_c_matrix1,1,min)*apply(IPCW_c_matrix2,1,min)

        I_K_ijv = signal[part]*weight_t[,part]*((cen_time_trt_last>(cen_time_con_last + tau_l))*event_con_last)
        g1Xc_K = apply(Xc_matrix,1,max)
        g2Xc_K = apply(cen_time_con_track,1,max)

        I_L_ijv = signal[part]*weight_c[,part]*((cen_time_con_last>(cen_time_trt_last + tau_l))*event_trt_last)
        g1Xt_L = apply(cen_time_trt_track,1,max)
        g2Xt_L = apply(Xt_matrix,1,max)

        status_KL[,(2*l-1)] = status_KL[,(2*l-1)] + get.Kijv(num1,num2,trt_id, con_id, ind_trt, ind_con, IPCW_t, I_K_ijv, g1Xc_K, g2Xc_K,
                                                             cen_time_trt, cen_time_con, cen_time_trt_long, cen_time_con_long,
                                                             event_trt, event_con, event_trt_long, event_con_long,
                                                             compare.vec_trt, compare.vec_con)
        status_KL[,(2*l)] = status_KL[,(2*l)] + get.Lijv(num1,num2,trt_id, con_id, ind_trt, ind_con, IPCW_c, I_L_ijv, g1Xt_L, g2Xt_L,
                                                         cen_time_trt, cen_time_con, cen_time_trt_long, cen_time_con_long,
                                                         event_trt, event_con, event_trt_long, event_con_long,
                                                         compare.vec_trt, compare.vec_con)
      }
    }

    Tie <- NA
    if(l == n_ep && l>1){
      Tie_t = -status_KL[,(2*l-1)]
      Tie_c = -status_KL[,(2*l)]

      pm_ind_split = strsplit(pm_ind,split = "")
      for(part in 1:(2^(l-1))){
        part_pm = pm_ind_split[[part]]
        IPCW_t_matrix1 = csurv1.x_con_mtau[con_id,l]; IPCW_c_matrix2 = csurv2.x_trt_mtau[trt_id,l]
        IPCW_t_matrix2 = csurv2.x_con[con_id,]; IPCW_c_matrix1 = csurv1.x_trt[trt_id,]
        Xt_matrix = cen_time_trt_track[,l] - tau_l
        Xc_matrix = cen_time_con_track[,l] - tau_l
        for(pmi in 1:(l-1)){
          if(part_pm[pmi] == "m"){
            IPCW_t_matrix1 = cbind(IPCW_t_matrix1,csurv1.x_con_mtau[con_id,l-pmi])
            IPCW_c_matrix2 = cbind(IPCW_c_matrix2,csurv2.x_trt_mtau[trt_id,l-pmi])
            Xt_matrix = cbind(Xt_matrix,cen_time_trt_track[,l-pmi] - tau_l)
            Xc_matrix = cbind(Xc_matrix,cen_time_con_track[,l-pmi] - tau_l)
          }else{
            IPCW_t_matrix1 = cbind(IPCW_t_matrix1,csurv1.x_con_ptau[con_id,l-pmi])
            IPCW_c_matrix2 = cbind(IPCW_c_matrix2,csurv2.x_trt_ptau[trt_id,l-pmi])
            Xt_matrix = cbind(Xt_matrix,cen_time_trt_track[,l-pmi] + tau_l)
            Xc_matrix = cbind(Xc_matrix,cen_time_con_track[,l-pmi] + tau_l)
          }
        }
        IPCW_t = apply(IPCW_t_matrix1,1,min)*apply(IPCW_t_matrix2,1,min)
        IPCW_c = apply(IPCW_c_matrix1,1,min)*apply(IPCW_c_matrix2,1,min)

        I_K_ijv = signal[part]*weight_t[,part]*((cen_time_trt_last>(cen_time_con_last - tau_l))*event_con_last)
        g1Xc_K = apply(Xc_matrix,1,max)
        g2Xc_K = apply(cen_time_con_track,1,max)

        I_L_ijv = signal[part]*weight_c[,part]*((cen_time_con_last>(cen_time_trt_last - tau_l))*event_trt_last)
        g1Xt_L = apply(cen_time_trt_track,1,max)
        g2Xt_L = apply(Xt_matrix,1,max)

        Tie_t = Tie_t + get.Kijv(num1,num2,trt_id, con_id, ind_trt, ind_con, IPCW_t, I_K_ijv, g1Xc_K, g2Xc_K,
                                 cen_time_trt, cen_time_con, cen_time_trt_long, cen_time_con_long,
                                 event_trt, event_con, event_trt_long, event_con_long,
                                 compare.vec_trt, compare.vec_con)
        Tie_c = Tie_c + get.Lijv(num1,num2,trt_id, con_id, ind_trt, ind_con, IPCW_c, I_L_ijv, g1Xt_L, g2Xt_L,
                                 cen_time_trt, cen_time_con, cen_time_trt_long, cen_time_con_long,
                                 event_trt, event_con, event_trt_long, event_con_long,
                                 compare.vec_trt, compare.vec_con)
        Tie = (Tie_t+Tie_c)/2
      }
    }
  }

  K = apply(as.matrix(status_KL[,seq(1,(2*n_ep-1),2)],ncol = n_ep),1,sum)
  L = apply(as.matrix(status_KL[,seq(2,(2*n_ep),2)],ncol = n_ep),1,sum)

  KL = cbind(trt_con$stratum, trt_con$pid_trt, trt_con$pid_con, K, L)
  colnames(KL) = c("stratum","pid_trt","pid_con","K","L")
  colnames(status_KL) = colnames(win_status)
  return(list(KL = as.data.frame(KL),status_KL = status_KL,Tie = Tie))
}
