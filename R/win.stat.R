win.stat<-function(data, ep_type, Z_t_trt = NULL, Z_t_con = NULL, iptw.weight = NULL, arm.name = c(1,2),
                   priority = c(1,2), alpha = 0.05, digit = 5, tau = 0, np_direction = "larger", horizon= Inf,
                   win.strategy = NULL,
                   pvalue = c("one-sided","two-sided"),
                   stratum.weight = c("unstratified","MH-type","wt.stratum1","wt.stratum2","equal"),
                   method = c("unadjusted","ipcw_tau","ipcw","covipcw","iptw","iptw_ipcw","iptw_covipcw"),
                   summary.print = TRUE, ...){
  #### match the argument
  pvalue = match.arg(pvalue)
  stratum.weight = match.arg(stratum.weight)
  method = match.arg(method)

  #### Remove missing values
  colname.ds = colnames(data)

  if(sum(is.na(data))>0){
    if(max(c("arm","trt","treat","treatment")%in%colname.ds)==TRUE){
      arm0 = data[,which(colname.ds%in%c("arm","trt","treat","treatment"))]
    }else{
      stop("The treatment variable is not found. Please rename the treatment variable to arm, trt, treat or treatment.")
    }
    ind.missing.trt = which(apply(data[arm0==arm.name[1],], 1, func<-function(x) sum(is.na(x))>0))
    ind.missing.con = which(apply(data[arm0==arm.name[2],], 1, func<-function(x) sum(is.na(x))>0))
    if(is.null(Z_t_trt) == FALSE){
      if("id"%in%colname.ds==TRUE){
        if(length(ind.missing.trt) > 0){
          id_trt0 = data[arm0==arm.name[1],which(colname.ds%in%c("id"))]
          Z_t_trt = Z_t_trt[Z_t_trt$id %in% id_trt0[-ind.missing.trt],]
        }
        if(length(ind.missing.con) > 0){
          id_con0 = data[arm0==arm.name[2],which(colname.ds%in%c("id"))]
          Z_t_con = Z_t_con[Z_t_con$id %in% id_con0[-ind.missing.con],]
        }
      }else{
        stop("The id variable is not found in Z_t_trt and Z_t_con.")
      }
    }
    data = na.omit(data)
    cat(length(ind.missing.trt)," and ",length(ind.missing.con),
        " objects with missing values are removed in the treatment and control group, respectively.","\n")
    warning("All the data entered of these objects will be removed due to incomplete/missing data.")
  }

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

  #### If np_direction is input as scalar, treat all the np_direction as the same.
  if(length(np_direction)==1){
    np_direction = rep(np_direction,n_ep)
  }

  if("smaller"%in%np_direction & method %in% c("ipcw","covipcw","ipcw_tau","iptw_ipcw","iptw_covipcw")){
    stop("The IPCW-adjusted approach is not applicable if smaller is specified for any endpoint in np_direction. Please try another method.")
  }

  #############################################################################################
  #### Reorganize the data
  #############################################################################################
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

  if(("stratum"%in%colname.ds)==TRUE && stratum.weight != "unstratified"){
    stratum = data[,which(colname.ds=="stratum")]
  }else{
    # The unstratified win statistics are calculated as the special case for stratified analysis with only one stratum.
    stratum = rep(1,n_total)
    cat("This analysis is unstratified.","\n")
  }

  Y = as.matrix(data[,which(stringr::str_detect(colname.ds,"Y"))])

  Delta = matrix(1,n_total,n_ep)
  if(max(c(stringr::str_detect(colname.ds,"Delta"),stringr::str_detect(colname.ds,"delta")))>0){
    ind_delta = which(stringr::str_detect(colname.ds,"Delta")|stringr::str_detect(colname.ds,"delta"))
    Delta[,which(ep_type=="tte")] = as.matrix(data[,ind_delta])
  }else if("tte"%in%ep_type){
    warning("No event status information detected. The default value of 1 is assigned to all individuals.")
  }

  Y[,which(ep_type%in%c("tte","continuous"))] = pmin(Y[,which(ep_type%in%c("tte","continuous"))], horizon)
  Delta[,which(ep_type%in%c("tte","continuous"))] = ifelse(Y[,which(ep_type%in%c("tte","continuous"))]==horizon,
                                                           1,
                                                           Delta[,which(ep_type%in%c("tte","continuous"))])

  df = data.frame(arm = arm,stratum = stratum,Delta = Delta,Y = Y)
  colnames(df) = c("arm","stratum",paste0("Delta_",1:n_ep),paste0("Y_",1:n_ep))
  rm(arm,stratum,Delta,Y)


  #############################################################################################
  #### Calculate the win statistics and make inference
  #############################################################################################
  res = switch (method,
                "unadjusted" = unadjusted.win.stat(df = df,n_total = n_total,arm.name = arm.name,id_trt=id_trt,id_con=id_con,
                                                   ep_type = ep_type,n_ep = n_ep,
                                                   priority = priority,tau = tau,np_direction = np_direction,
                                                   win.strategy = win.strategy,
                                                   alpha = alpha,digit = digit,pvalue = pvalue,stratum.weight = stratum.weight,
                                                   summary.print = summary.print, ...),
                "ipcw_tau" = ipcw.adjusted.tau.win.stat(df = df,n_total = n_total,arm.name = arm.name,id_trt=id_trt,id_con=id_con,
                                                        ep_type = ep_type,n_ep = n_ep,priority = priority,tau = tau,
                                                        np_direction = np_direction,win.strategy = win.strategy,
                                                        alpha = alpha,digit = digit,pvalue = pvalue,stratum.weight = stratum.weight,
                                                        summary.print = summary.print, ...),
                "ipcw" = ipcw.win.stat(df = df,n_total = n_total,arm.name = arm.name,id_trt=id_trt,id_con=id_con,
                                       ep_type = ep_type,n_ep = n_ep,
                                       priority = priority,tau = tau,np_direction = np_direction,win.strategy = win.strategy,
                                       alpha = alpha,digit = digit,pvalue = pvalue,stratum.weight = stratum.weight,
                                       summary.print = summary.print, ...),
                "covipcw" = covipcw.win.stat(df = df,Z_t_trt = Z_t_trt,Z_t_con = Z_t_con,n_total = n_total,
                                             arm.name = arm.name,id_trt=id_trt,id_con=id_con,ep_type = ep_type,n_ep = n_ep,
                                             priority = priority,tau = tau,np_direction = np_direction,win.strategy = win.strategy,
                                             alpha = alpha,digit = digit,pvalue = pvalue,stratum.weight = stratum.weight,
                                             summary.print = summary.print, ...),
                "iptw" = iptw.adjusted.win.stat(df = df,iptw.weight = iptw.weight,n_total = n_total,
                                                arm.name = arm.name,id_trt=id_trt,id_con=id_con,ep_type = ep_type,n_ep = n_ep,
                                                priority = priority,tau = tau,np_direction = np_direction,win.strategy = win.strategy,
                                                alpha = alpha,digit = digit,pvalue = pvalue,stratum.weight = stratum.weight,
                                                summary.print = summary.print, ...),
                "iptw_ipcw" = iptw_ipcw.win.stat(df = df,iptw.weight = iptw.weight,n_total = n_total,arm.name = arm.name,
                                                 id_trt=id_trt,id_con=id_con,ep_type = ep_type,n_ep = n_ep,priority = priority,
                                                 tau = tau,np_direction = np_direction,win.strategy = win.strategy,alpha = alpha,
                                                 digit = digit,pvalue = pvalue,stratum.weight = stratum.weight,
                                                 summary.print = summary.print, ...),
                "iptw_covipcw" = iptw_covipcw.win.stat(df = df,Z_t_trt = Z_t_trt,Z_t_con = Z_t_con,iptw.weight = iptw.weight,
                                                       n_total = n_total,arm.name = arm.name,id_trt=id_trt,id_con=id_con,
                                                       ep_type = ep_type,n_ep = n_ep,priority = priority,tau = tau,
                                                       np_direction = np_direction,win.strategy = win.strategy,
                                                       alpha = alpha,digit = digit,pvalue = pvalue,stratum.weight = stratum.weight,
                                                       summary.print = summary.print, ...)
  )

  return(res)
}
