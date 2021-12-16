sim.data<-function(randomseed = 12345, n_trt = 200, n_con = 200, n_ep = 2, n_stratum = 1,
                   arm.name = c(1,2), ep_type, cdist.rate, sim_method = "copula",
                   copula_trt = NULL, margins_trt = NULL, paramMargins_trt = NULL,
                   copula_con = NULL, margins_con = NULL, paramMargins_con = NULL,
                   rate_trt = NULL, rate_con = NULL, max_accrual_time = NULL){
  #### set random seed
  set.seed(randomseed)

  #### If ep_type is input as scalar, treat all the endpoints as the same type.
  if(length(ep_type)==1){
    print(paste0("The outcome type for all the endpoints: ",ep_type))
    ep_type = rep(ep_type,n_ep)
  }
  n_tte = sum(ep_type=="tte")

  #### Generate the vector for arm and stratum
  arm = c(rep(arm.name[1],n_trt),rep(arm.name[2],n_con))
  stratum = sample(1:n_stratum,(n_trt+n_con),replace = TRUE)


  #### Generate the outcome
  if(sim_method == "copula"){
    my.mvd_trt = copula::mvdc(copula = copula_trt, margins = margins_trt, paramMargins = paramMargins_trt)
    Y_trt = copula::rMvdc(n_trt,my.mvd_trt)

    my.mvd_con = copula::mvdc(copula = copula_con, margins = margins_con, paramMargins = paramMargins_con)
    Y_con = copula::rMvdc(n_con,my.mvd_con)

    Y = rbind(Y_trt,Y_con)

    if(n_tte>0){
      if(length(cdist.rate)==1){
        cdist.rate = rep(cdist.rate,n_tte)
      }

      if(length(cdist.rate)!=n_tte){
        stop("The length of cdist.rate does not match the number of TTE endpoints.")
      }

      ind.tte = which(ep_type=="tte")
      Delta = NULL
      for(ci in 1:n_tte){
        Ctime = rexp(n = (n_trt+n_con), rate = cdist.rate[ci])
        Delta = cbind(Delta, 1 * (Ctime > Y[,ind.tte[ci]]))
        Y_temp = Y[,ind.tte[ci]] * (Ctime > Y[,ind.tte[ci]]) + Ctime * (Ctime <= Y[,ind.tte[ci]])
        Y[,ind.tte[ci]] = Y_temp
      }

      sim.data = data.frame(arm = arm,stratum = stratum,Delta = Delta,Y = Y)
      colnames(sim.data) = c("arm","stratum",paste0("Delta_",ind.tte),paste0("Y_",1:n_ep))
    }else{
      sim.data = data.frame(arm = arm,stratum = stratum,Y = Y)
      colnames(sim.data) = c("arm","stratum",paste0("Y_",1:n_ep))
    }
  }else if(sim_method == "tte_exponential"){
    #### Return error if not designed endpoints (two tte endpoints for progression and death time)
    if(max(ep_type!="tte")==TRUE | n_ep != 2){
      stop("The option tte_exponential is designed only for two TTE endpoints with the more important TTE endpoint expected to occur later.")
    }
    if(length(cdist.rate)==1){
      cdist.rate = rep(cdist.rate,n_tte)
    }

    Ctime1_trt = rexp(n = n_trt, rate = cdist.rate[1])
    Ctime2_trt = rexp(n = n_trt, rate = cdist.rate[2])

    time1_trt = rexp(n_trt, rate = rate_trt[1])
    time2_trt = rexp(n_trt, rate = rate_trt[2])

    endpoint1_trt = apply(cbind(time1_trt,time2_trt,Ctime1_trt,Ctime2_trt), 1, min)
    endpoint2_trt = apply(cbind(time2_trt,Ctime2_trt), 1, min)

    Delta_1_trt = 1*(apply(cbind(time2_trt, Ctime1_trt,Ctime2_trt),1,min) > time1_trt)
    Delta_2_trt = 1*(Ctime2_trt > time2_trt)

    Ctime1_con = rexp(n = n_con, rate = cdist.rate[1])
    Ctime2_con = rexp(n = n_con, rate = cdist.rate[2])

    time1_con = rexp(n_con, rate = rate_con[1])
    time2_con = rexp(n_con, rate = rate_con[2])

    endpoint1_con = apply(cbind(time1_con,time2_con,Ctime1_con,Ctime2_con), 1, min)
    endpoint2_con = apply(cbind(time2_con,Ctime2_con), 1, min)

    Delta_1_con = 1*(apply(cbind(time2_con, Ctime1_con,Ctime2_con),1,min) > time1_con)
    Delta_2_con = 1*(Ctime2_con > time2_con)

    Y = rbind(cbind(endpoint1_trt,endpoint2_trt),cbind(endpoint1_con,endpoint2_con))

    Delta = rbind(cbind(Delta_1_trt,Delta_2_trt),cbind(Delta_1_con,Delta_2_con))

    sim.data = data.frame(arm = arm,stratum = stratum,Delta = Delta,Y = Y)
    colnames(sim.data) = c("arm","stratum",paste0("Delta_",c(1,2)),paste0("Y_",c(1,2)))
  }

  if(!is.null(max_accrual_time)){
    Start_time = runif((n_trt+n_con), min = 0, max = max_accrual_time)
    sim.data = cbind(sim.data, Start_time = Start_time)
  }

  return(sim.data)
}
