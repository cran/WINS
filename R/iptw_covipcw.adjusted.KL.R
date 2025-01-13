iptw_covipcw.adjusted.KL<-function(win_status, trt, con, trt_con, Z_t_trt, Z_t_con, priority, n_ep, ep_type, weight.trt, weight.con){
  #### Obtain the number of treatment and control
  num1 = dim(trt)[1]
  num2 = dim(con)[1]

  #### Calculate K and L
  status_KL = win_status

  trt_id = as.factor(trt_con$pid_trt); levels(trt_id) = 1:num1
  con_id = as.factor(trt_con$pid_con); levels(con_id) = 1:num2

  weight.KL = weight.trt[trt_id]*weight.con[con_id]

  status_KL = win_status*weight.KL

  for(l in n_ep){
    ind.prior = priority[l]
    if(ep_type[ind.prior]=="tte"){
      # Setup survival object
      cen_time_trt = trt[,(2+n_ep+ind.prior)]
      cen_time_con = con[,(2+n_ep+ind.prior)]

      event_trt = trt[,(2+ind.prior)]
      event_con = con[,(2+ind.prior)]

      ori_trt = data.frame(id = trt$pid_trt,ctime = cen_time_trt)
      ori_con = data.frame(id = con$pid_con,ctime = cen_time_con)

      long_trt = survival::tmerge(ori_trt,ori_trt,id=id,tstop = ctime)
      long_trt = survival::tmerge(long_trt,Z_t_trt,id=id,event = event(time))
      long_trt$event = (1-long_trt$event)*(event_trt[long_trt$id]==0)
      long_trt = as.data.frame(cbind(long_trt[,-1],Z_t_trt[,-c(1,2)]))
      colnames(long_trt) = c("ctime","tstart","tstop","event",paste0("Z",1:(ncol(Z_t_trt)-2)))

      long_con = survival::tmerge(ori_con,ori_con,id=id,tstop = ctime)
      long_con = survival::tmerge(long_con,Z_t_con,id=id,event = event(time))
      long_con$event = (1-long_con$event)*(event_con[long_con$id]==0)
      long_con = as.data.frame(cbind(long_con[,-1],Z_t_con[,-c(1,2)]))
      colnames(long_con) = c("ctime","tstart","tstop","event",paste0("Z",1:(ncol(Z_t_con)-2)))

      # Obtain the estimator of the survival function for censoring time using time-dependent Cox model
      if(!ncol(Z_t_trt)==ncol(Z_t_con)){
        stop("The number of covariates is different in treatment and control group!")
      }

      csurv1.fit = survival::survfit(survival::coxph(survival::Surv(tstart, tstop, event)~
                                                       as.matrix(long_trt[,-c(1:4)]),
                                                     data = long_trt),weights = weight.trt[long_trt$id])
      csurv1.fit$time = c(0,csurv1.fit$time)
      csurv1.fit$surv = c(1,csurv1.fit$surv)
      L.tpt_trt = length(csurv1.fit$time)

      csurv1.indx_trt = apply(cen_time_trt >= t(array(rep(csurv1.fit$time, num1), c(L.tpt_trt, num1))), 1, sum)
      csurv1.x_trt = csurv1.fit$surv[csurv1.indx_trt]*as.numeric((event_trt==1))+(1-as.numeric((event_trt==1)))

      csurv1.indx_con = apply(cen_time_con >= t(array(rep(csurv1.fit$time, num2), c(L.tpt_trt, num2))), 1, sum)
      csurv1.x_con = csurv1.fit$surv[csurv1.indx_con]*as.numeric((event_con==1))+(1-as.numeric((event_con==1)))

      csurv2.fit = survival::survfit(survival::coxph(survival::Surv(tstart, tstop, event)~
                                                       as.matrix(long_con[,-c(1:4)]),
                                                     data = long_con),weights = weight.con[long_con$id])
      csurv2.fit$time = c(0,csurv2.fit$time)
      csurv2.fit$surv = c(1,csurv2.fit$surv)
      L.tpt_con = length(csurv2.fit$time)

      csurv2.indx_trt = apply(cen_time_trt >= t(array(rep(csurv2.fit$time, num1), c(L.tpt_con, num1))), 1, sum)
      csurv2.x_trt = csurv2.fit$surv[csurv2.indx_trt]*as.numeric((event_trt==1))+(1-as.numeric((event_trt==1)))

      csurv2.indx_con = apply(cen_time_con >= t(array(rep(csurv2.fit$time, num2), c(L.tpt_con, num2))), 1, sum)
      csurv2.x_con = csurv2.fit$surv[csurv2.indx_con]*as.numeric((event_con==1))+(1-as.numeric((event_con==1)))

      trt_id = as.factor(trt_con$pid_trt); levels(trt_id) = 1:num1
      con_id = as.factor(trt_con$pid_con); levels(con_id) = 1:num2

      IPCW1 = csurv1.x_trt[trt_id]*csurv2.x_trt[trt_id]
      IPCW2 = csurv1.x_con[con_id]*csurv2.x_con[con_id]

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
