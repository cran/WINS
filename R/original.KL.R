original.KL<-function(win_status, trt_con, n_ep){
  K = apply(as.matrix(win_status[,seq(1,(2*n_ep-1),2)],ncol = n_ep),1,sum)
  L = apply(as.matrix(win_status[,seq(2,(2*n_ep),2)],ncol = n_ep),1,sum)

  KL = cbind(trt_con$stratum, trt_con$pid_trt, trt_con$pid_con, K, L)
  colnames(KL) = c("stratum","pid_trt","pid_con","K","L")
  return(list(KL = as.data.frame(KL),status_KL = win_status))
}
