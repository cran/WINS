get.Lijv <- function(num1,num2,trt_id, con_id, ind_trt, ind_con, IPCW_c, I_L_ijv, g1Xt_L, g2Xt_L,
                     cen_time_trt, cen_time_con, cen_time_trt_long, cen_time_con_long,
                     event_trt, event_con, event_trt_long, event_con_long,
                     compare.vec_trt, compare.vec_con){
  Lijv_1 = I_L_ijv/IPCW_c
  temp_sum_1 = NULL
  for(i in 1:num1){
    temp_sum_1 = c(temp_sum_1,sum(I_L_ijv[trt_id == i]))
  }
  temp_sum_1_long = temp_sum_1[trt_id]
  temp_sum_2_long = compare.vec_con[con_id]
  Lijv_4 = ((event_con_long == 0)*(g2Xt_L >= cen_time_con_long)*temp_sum_1_long)/
    (IPCW_c * temp_sum_2_long)
  Lijv_5 = 0
  for(q in 1:num2){
    Lijv_5 = Lijv_5 +
      ((event_con[q] == 0)*(cen_time_con_long >= cen_time_con[q])*(g2Xt_L >= cen_time_con[q])*temp_sum_1_long)/
      (IPCW_c * (compare.vec_con[q]^2))
  }
  temp_sum_3_long = 0
  for(k in 1:num1){
    temp_sum_3_long = temp_sum_3_long +
      (I_L_ijv[c(((k-1)*num2+1):(k*num2))][con_id])*
      (g1Xt_L[(k-1)*num2+1] >= cen_time_trt_long)/IPCW_c
  }
  temp_sum_4_long = compare.vec_trt[trt_id]
  Lijv_2 = ((event_trt_long == 0)*temp_sum_3_long)/temp_sum_4_long
  Lijv_3 = 0
  for(p in 1:num1){
    Lijv_3 = Lijv_3 +
      ((event_trt[p] == 0)*(cen_time_trt_long >= cen_time_trt[p])*(temp_sum_3_long[c(((p-1)*num2+1):(p*num2))][con_id]))/
      (compare.vec_trt[p]^2)
  }

  Lijv = Lijv_1 + Lijv_2 - Lijv_3 + Lijv_4 - Lijv_5
  return(Lijv)
}
