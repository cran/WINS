iptw.adjusted.KL<-function(win_status, trt, con, trt_con, priority, n_ep, ep_type, weight.trt, weight.con){
  #### Obtain the number of treatment and control
  num1 = dim(trt)[1]
  num2 = dim(con)[1]

  #### Calculate K and L
  status_KL = win_status

  trt_id = as.factor(trt_con$pid_trt); levels(trt_id) = 1:num1
  con_id = as.factor(trt_con$pid_con); levels(con_id) = 1:num2

  weight.KL = weight.trt[trt_id]*weight.con[con_id]

  status_KL = win_status*weight.KL

  K = apply(as.matrix(status_KL[,seq(1,(2*n_ep-1),2)],ncol = n_ep),1,sum)
  L = apply(as.matrix(status_KL[,seq(2,(2*n_ep),2)],ncol = n_ep),1,sum)

  KL = cbind(trt_con$stratum, trt_con$pid_trt, trt_con$pid_con, K, L)
  colnames(KL) = c("stratum","pid_trt","pid_con","K","L")
  return(list(KL = as.data.frame(KL),status_KL = status_KL))
}
