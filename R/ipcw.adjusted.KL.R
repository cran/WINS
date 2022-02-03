ipcw.adjusted.KL<-function(win_status, trt, con, trt_con, priority, n_ep, ep_type){
  #### Obtain the number of treatment and control
  num1 = dim(trt)[1]
  num2 = dim(con)[1]

  #### Calculate K and L
  status_KL = win_status

  for(l in n_ep){
    ind.prior = priority[l]
    if(ep_type[ind.prior]=="tte"){
      # Setup survival object
      cen_time_trt = trt[,(2+n_ep+ind.prior)]
      cen_time_con = con[,(2+n_ep+ind.prior)]
      ##cen_time = c(cen_time_trt,cen_time_con)

      event_trt = trt[,(2+ind.prior)]
      event_con = con[,(2+ind.prior)]
      ##event = c(event_trt,event_con)

      # Obtain Kaplan-Meier estimator, the survival rate is truncated by 0.1 if small in tail.
      surv_trt = survival::Surv(time = cen_time_trt, event = as.numeric(event_trt==0))
      csurv1.fit = survival::survfit(surv_trt~1,type='kaplan-meier')
      ind.trun1 = (csurv1.fit$surv >= 0.1)
      csurv1.fit$time = c(0,csurv1.fit$time[ind.trun1])
      csurv1.fit$surv = c(1,csurv1.fit$surv[ind.trun1])
      L.tpt_trt = length(csurv1.fit$time)

      csurv1.indx_trt = apply(cen_time_trt >= t(array(rep(csurv1.fit$time, num1), c(L.tpt_trt, num1))), 1, sum)
      csurv1.x_trt = csurv1.fit$surv[csurv1.indx_trt]*as.numeric((event_trt==1))+(1-as.numeric((event_trt==1)))

      csurv1.indx_con = apply(cen_time_con >= t(array(rep(csurv1.fit$time, num2), c(L.tpt_trt, num2))), 1, sum)
      csurv1.x_con = csurv1.fit$surv[csurv1.indx_con]*as.numeric((event_con==1))+(1-as.numeric((event_con==1)))

      surv_con = survival::Surv(time = cen_time_con, event = as.numeric(event_con==0))
      csurv2.fit = survival::survfit(surv_con~1,type='kaplan-meier')
      ind.trun2 = (csurv2.fit$surv >= 0.1)
      csurv2.fit$time = c(0,csurv2.fit$time[ind.trun2])
      csurv2.fit$surv = c(1,csurv2.fit$surv[ind.trun2])
      L.tpt_con = length(csurv2.fit$time)

      csurv2.indx_trt = apply(cen_time_trt >= t(array(rep(csurv2.fit$time, num1), c(L.tpt_con, num1))), 1, sum)
      csurv2.x_trt = csurv2.fit$surv[csurv2.indx_trt]*as.numeric((event_trt==1))+(1-as.numeric((event_trt==1)))

      csurv2.indx_con = apply(cen_time_con >= t(array(rep(csurv2.fit$time, num2), c(L.tpt_con, num2))), 1, sum)
      csurv2.x_con = csurv2.fit$surv[csurv2.indx_con]*as.numeric((event_con==1))+(1-as.numeric((event_con==1)))

      IPCW1 = csurv1.x_trt[trt_con$pid_trt]*csurv2.x_trt[trt_con$pid_trt]
      IPCW2 = csurv1.x_con[(trt_con$pid_con)]*csurv2.x_con[(trt_con$pid_con)]

      status_KL[,(2*l-1)] = win_status[,(2*l-1)]/IPCW2
      status_KL[,(2*l)] = win_status[,(2*l)]/IPCW1
    }
  }

  K = apply(as.matrix(status_KL[,seq(1,(2*n_ep-1),2)],ncol = n_ep),1,sum)
  L = apply(as.matrix(status_KL[,seq(2,(2*n_ep),2)],ncol = n_ep),1,sum)

  KL = cbind(trt_con$stratum, trt_con$pid_trt, trt_con$pid_con, K, L)
  colnames(KL) = c("stratum","pid_trt","pid_con","K","L")
  return(list(KL = as.data.frame(KL),status_KL = status_KL))
}
