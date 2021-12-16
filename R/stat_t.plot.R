stat_t.plot<-function(data, Ctime = Inf, arm.name = c(1,2), priority = c(1,2),
                        statistic = c("WR","NB","WO"), Z_t_trt = NULL, Z_t_con = NULL,tau = 0,
                        weight = c("unstratified","MH-type","wt.stratum1","wt.stratum2","equal"),
                        censoring_adjust = c("No","IPCW","CovIPCW"),
                        win.strategy = NULL, plotTimeUnit = NULL, plot_CI = FALSE, alpha = 0.05,...){
  #### match the argument
  statistic = match.arg(statistic)
  weight = match.arg(weight)
  censoring_adjust = match.arg(censoring_adjust)

  #### obtain the number of endpoints and total number of individuals
  n_ep = length(priority)
  n_total = dim(data)[1]

  ep_type = rep("tte",n_ep)

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

  if(("stratum"%in%colname.ds)==TRUE && weight != "unstratified"){
    stratum = data[,which(colname.ds=="stratum")]
  }else{
    # The unstratified win statistics are calculated as the special case for stratified analysis with only one stratum.
    stratum = rep(1,n_total)
    cat("This analysis is unstratified.","\n")
  }

  if(max(c("Start_time","start_time","entry_time","Entry_time")%in%colname.ds)==TRUE){
    Start_time = data[,which(colname.ds%in%c("Start_time","start_time","entry_time","Entry_time"))]
  }else{
    Start_time = rep(0,nrow(data))
    warning("The study entry time is missing, by default zero will be assigned to all subjects.")
  }

  Delta = matrix(1,n_total,n_ep)
  if(max(stringr::str_detect(colname.ds,"Delta"))>0){
    Delta[,which(ep_type=="tte")] = as.matrix(data[,which(stringr::str_detect(colname.ds,"Delta"))])
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

  #############################################################################################
  #### print the plot
  #############################################################################################
  #### obtain the win statistic for each time point
  if(!plot_CI){
    win_stat_t = NULL
    for(j in 1:length(Ctime)){
      res_t = get.win.stat_t(trt = trt, con = con, ep_type = ep_type, priority = priority,
                             Ctimej = Ctime[j], Start_time_trt = Start_time_trt,
                             Start_time_con = Start_time_con, Z_t_trt = Z_t_trt, Z_t_con = Z_t_con,
                             tau = tau, weight = weight, censoring_adjust = censoring_adjust,
                             win.strategy = win.strategy, ...)
      ind_stat = switch (statistic,
                         "WR" = 1,
                         "NB" = 2,
                         "WO" = 3
      )
      win_stat_t = c(win_stat_t, res_t[[ind_stat]])
    }
    df1 = data.frame(time = Ctime, win_stat = win_stat_t)
  }else{
    win_stat_t = NULL
    for(j in 1:length(Ctime)){
      res_t = get.win.stat_t(trt = trt, con = con, ep_type = ep_type, priority = priority,
                             Ctimej = Ctime[j], Start_time_trt = Start_time_trt,
                             Start_time_con = Start_time_con, Z_t_trt = Z_t_trt, Z_t_con = Z_t_con,
                             tau = tau, weight = weight, censoring_adjust = censoring_adjust,
                             win.strategy = win.strategy, return_CI = TRUE, pvalue = pvalue,
                             alpha = alpha, ...)
      ind_stat = switch (statistic,
                         "WR" = 1,
                         "NB" = 2,
                         "WO" = 3
      )
      win_stat_t = rbind(win_stat_t, c(res_t$Win_statisitc[[ind_stat]], res_t$CI_t[[ind_stat]]))
    }
    colnames(win_stat_t) = c("win_stat", "lower_ci", "upper_ci")

    df1 = data.frame(time = Ctime, win_stat_t)
  }

  if(!is.finite(max(Ctime))){
    warning("With Ctime set as infinite, no plot will be output. The usual win statistics are returned.")
    return(list(statistic = statistic, time = Ctime, win_stat = win_stat_t))
  }else{
    #### plot the win statistic over time
    x_lab_name = ifelse(is.null(plotTimeUnit),
                        yes = "Study time",
                        no = paste0("Study time in ", plotTimeUnit))
    y_lab_name = switch (statistic,
                         "WR" = "Win Ratio",
                         "NB" = "Net Benefit",
                         "WO" = "Win Odds"
    )
    g1 = ggplot2::ggplot(data = subset(df1, time<=max_study_time), ggplot2::aes(x=time, y=win_stat, group=1)) +
      ggplot2::geom_line(color="blue") +
      ggplot2::geom_point(color="orange") +
      ggplot2::xlab(x_lab_name) +
      ggplot2::ylab(y_lab_name) +
      ggplot2::xlim(min(Ctime),min(max_study_time+1,max(Ctime))) +
      ggplot2::ylim(min(c(0, min(df1[,-1]))),1.1*max(df1[,-1])) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(color = "grey20", size = 8, angle = 0, hjust = .5, vjust = .5, face = "plain"),
            axis.text.y = ggplot2::element_text(color = "grey20", size = 8, angle = 0, hjust = 1, vjust = 0, face = "plain"),
            axis.title.x = ggplot2::element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
            axis.title.y = ggplot2::element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"))

    if(max(Ctime)>=max_study_time){
      g1 = g1 + ggplot2::geom_vline(xintercept = max_study_time, linetype="dotted", color = "grey", size=1.2) +
        ggplot2::annotate(geom="text", x=0.8*max_study_time, y= 1.1*max(df1[,-1]), label="Maximum study time", color="grey")
    }

    if(plot_CI){
      g1 = g1 + ggplot2::geom_ribbon(ggplot2::aes(ymin=lower_ci, ymax=upper_ci), alpha=0.2)
    }

    print(g1)

    values = data.frame(time = Ctime[1:ind.max],win_stat_t[1:ind.max,])

    return(list(statistic = statistic, values = values))
  }
}

