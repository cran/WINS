partition_t.plot<-function(data, Ctime = Inf, arm.name = c(1,2), priority = c(1,2), censoring_adjust = "No",
                           Z_t_trt = NULL, Z_t_con = NULL,tau = 0, np_direction = "larger",
                           plotTimeUnit = NULL, trt_group = c("both","trt","con"),
                           win.strategy = NULL, ...){
  #### match the argument
  trt_group = match.arg(trt_group)

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
  }

  #### obtain the number of endpoints and total number of individuals
  n_ep = length(priority)
  n_total = dim(data)[1]

  ep_type = rep("tte",n_ep)

  #### If tau is input as scalar, treat all the tau as the same.
  if(length(tau)==1){
    tau = rep(tau,n_ep)
  }

  #### If np_direction is input as scalar, treat all the np_direction as the same.
  if(length(np_direction)==1){
    np_direction = rep(np_direction,n_ep)
  }

  #############################################################################################
  #### Reorganize the data
  #############################################################################################
  if(max(c("arm","trt","treat","treatment")%in%colname.ds)==TRUE){
    arm = data[,which(colname.ds%in%c("arm","trt","treat","treatment"))]
  }else{
    stop("The treatment variable is not found. Please rename the treatment variable to arm, trt, treat or treatment.")
  }

  if(("stratum"%in%colname.ds)==TRUE){
    stratum = data[,which(colname.ds=="stratum")]
  }else{
    stratum = rep(1,n_total)
  }

  if(max(c("Start_time","start_time","entry_time","Entry_time")%in%colname.ds)==TRUE){
    Start_time = data[,which(colname.ds%in%c("Start_time","start_time","entry_time","Entry_time"))]
  }else{
    Start_time = rep(0,nrow(data))
    warning("The study entry time is missing, by default zero will be assigned to all subjects.")
  }

  Delta = matrix(1,n_total,n_ep)
  if(max(c(stringr::str_detect(colname.ds,"Delta"),stringr::str_detect(colname.ds,"delta")))>0){
    ind_delta = which(stringr::str_detect(colname.ds,"Delta")|stringr::str_detect(colname.ds,"delta"))
    Delta[,which(ep_type=="tte")] = as.matrix(data[,ind_delta])
  }else if("tte"%in%ep_type){
    warning("No event status information detected. The default value of 1 is assigned to all individuals.")
  }

  Y = as.matrix(data[,which(stringr::str_detect(colname.ds,"Y"))])

  #### obtain the maximum study time
  study_time = sweep(Y[,which(ep_type=="tte")], 1, Start_time, "+")
  max_study_time = max(study_time)
  ind.max = max(which(Ctime<=max_study_time))

  data_mix = data.frame(arm = arm,stratum = stratum,Delta = Delta,Y = Y)
  colnames(data_mix) = c("arm","stratum",paste0("Delta_",1:n_ep),paste0("Y_",1:n_ep))
  rm(arm,stratum,Delta,Y)

  #### obtain the number of treatment and control
  n_trt = sum(data_mix$arm==arm.name[1]); n_con = sum(data_mix$arm==arm.name[2])

  #### obtain the start time for each individual in the treatment and control group
  Start_time_trt = Start_time[data_mix$arm==arm.name[1]]
  Start_time_con = Start_time[data_mix$arm==arm.name[2]]

  #############################################################################################
  #### obtain the unmatched pairs
  #############################################################################################
  #### pair the individuals in the treatment and control group
  trt = data.frame(1:n_trt,data_mix[data_mix$arm==arm.name[1],-1])
  colnames(trt) = c("pid_trt","stratum",paste0(colnames(data_mix)[-c(1,2)],"_trt"))
  con = data.frame(1:n_con,data_mix[data_mix$arm==arm.name[2],-1])
  colnames(con) = c("pid_con","stratum",paste0(colnames(data_mix)[-c(1,2)],"_con"))

  #### obtain the win statistic for each time
  win_status_trt_t = NULL; win_status_con_t = NULL; win_prop_tie_t = NULL
  for(j in 1:length(Ctime)){
    res = get.win.stat_t(trt = trt, con = con, ep_type = ep_type, priority = priority,
                         Ctimej = Ctime[j], Start_time_trt = Start_time_trt, censoring_adjust = censoring_adjust,
                         Start_time_con = Start_time_con, Z_t_trt = Z_t_trt, Z_t_con = Z_t_con,
                         tau = tau, np_direction = np_direction, win.strategy = win.strategy, status_only = TRUE, ...)

    temp_trt = apply(res$win_status[,seq(1,2*n_ep-1,2)], 2, mean)
    prop_trt = cbind(time = rep(Ctime[j],n_ep), Endpoint = priority, proportion = temp_trt)

    win_status_trt_t = rbind(win_status_trt_t, prop_trt)

    temp_con = apply(res$win_status[,seq(2,2*n_ep,2)], 2, mean)
    prop_con = cbind(time = rep(Ctime[j],n_ep), Endpoint = priority, proportion = temp_con)

    win_status_con_t = rbind(win_status_con_t, prop_con)

    temp_tie = 1 - sum(temp_trt + temp_con)
    prop_tie = c(Ctime[j], temp_tie)

    win_prop_tie_t = rbind(win_prop_tie_t, prop_tie)
  }

  colnames(win_prop_tie_t) = c("time","proportion of ties")
  rownames(win_prop_tie_t) = 1:nrow(win_prop_tie_t)

  colnames(win_status_trt_t) = c("time","endpoint","proportion")
  win_status_trt_t = as.data.frame(win_status_trt_t)

  win_status_trt_t$endpoint <- factor(win_status_trt_t$endpoint , levels=c(1:n_ep))
  win_trt_t <- reshape2::dcast(win_status_trt_t, time~endpoint,value.var = "proportion")
  colnames(win_trt_t) =  c("time",paste0("endpoint",1:n_ep))

  colnames(win_status_con_t) = c("time","endpoint","proportion")
  win_status_con_t = as.data.frame(win_status_con_t)

  win_status_con_t$endpoint <- factor(win_status_con_t$endpoint , levels=c(1:n_ep))
  win_con_t <- reshape2::dcast(win_status_con_t, time~endpoint,value.var = "proportion")
  colnames(win_con_t) =  c("time",paste0("endpoint",1:n_ep))

  x_lab_name = ifelse(is.null(plotTimeUnit),
                      yes = "Study time",
                      no = paste0("Study time in ", plotTimeUnit))

  win_prop_ep_trt = win_status_trt_t$proportion[(nrow(win_status_trt_t)-(n_ep-1)):nrow(win_status_trt_t)]
  win_prop_ep_con = win_status_con_t$proportion[(nrow(win_status_con_t)-(n_ep-1)):nrow(win_status_con_t)]

  y_u = 1.1*max(c(sum(win_prop_ep_con),sum(win_prop_ep_trt)))
  y_prop_trt = cumsum(win_prop_ep_trt[order(priority,decreasing = TRUE)])
  y_prop_con = cumsum(win_prop_ep_con[order(priority,decreasing = TRUE)])

  g2 = ggplot2::ggplot(data = subset(win_status_trt_t, time<=max_study_time),
                       ggplot2::aes(x=time, y=proportion, fill=endpoint)) +
    ggplot2::geom_area(alpha=0.6 , size=.5, colour="white") +
    viridis::scale_fill_viridis(discrete = T) +
    ggplot2::xlab(x_lab_name) +
    ggplot2::ylab("Win Proportion") +
    ggplot2::xlim(min(win_status_trt_t$time),min(max_study_time+1,max(win_status_trt_t$time))) +
    ggplot2::ylim(0,y_u) +
    ggplot2::ggtitle("The win proportion for the treatment group") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(color = "grey20", size = 8, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = ggplot2::element_text(color = "grey20", size = 8, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = ggplot2::element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = ggplot2::element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          legend.title = ggplot2::element_text(color = "grey20", size = 6),
          legend.text = ggplot2::element_text(color = "grey20", size = 4))

  if(max(Ctime)>=max_study_time){
    g2 = g2 + ggplot2::geom_vline(xintercept = max_study_time, linetype="dotted", color = "grey", size=1.2) +
      ggplot2::annotate(geom="text", x=0.8*max_study_time, y= y_u, label="Maximum study time", color="grey")
  }

  if(trt_group == "both"){
    g1 = ggplot2::ggplot(data = subset(win_status_con_t, time<=max_study_time),
                         ggplot2::aes(x=time, y=proportion, fill=endpoint)) +
      ggplot2::geom_area(alpha=0.6 , size=.5, colour="white") +
      viridis::scale_fill_viridis(discrete = T) +
      ggplot2::xlab(x_lab_name) +
      ggplot2::ylab("Win Proportion") +
      ggplot2::xlim(min(win_status_con_t$time),min(max_study_time+1,max(win_status_con_t$time))) +
      ggplot2::ylim(0,y_u) +
      ggplot2::ggtitle("The win proportion for the control group") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(color = "grey20", size = 8, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = ggplot2::element_text(color = "grey20", size = 8, angle = 0, hjust = 1, vjust = 0, face = "plain"),
            axis.title.x = ggplot2::element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = ggplot2::element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            legend.title = ggplot2::element_text(color = "grey20", size = 6),
            legend.text = ggplot2::element_text(color = "grey20", size = 4))

    if(max(Ctime)>=max_study_time){
      g1 = g1 + ggplot2::geom_vline(xintercept = max_study_time, linetype="dotted", color = "grey", size=1.2) +
        ggplot2::annotate(geom="text", x=0.8*max_study_time, y=y_u, label="Maximum study time", color="grey")
    }

    print(ggpubr::ggarrange(g2, g1, ncol=2, common.legend = TRUE, legend="right"))
  }

  if(trt_group == "trt"){
    print(g2)
  }

  if(trt_group == "con"){
    print(g1)
  }

  return(list(win_trt_t = win_trt_t[1:ind.max,],
              win_con_t = win_con_t[1:ind.max,],
              win_tie_t = win_prop_tie_t[1:ind.max,],
              max_study_time = max_study_time))
}
