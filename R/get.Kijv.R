get.Kijv <- function(num1,num2,trt_id, con_id, ind_trt, ind_con, IPCW_t, I_K_ijv, g1Xc_K, g2Xc_K,
                     cen_time_trt, cen_time_con, cen_time_trt_long, cen_time_con_long,
                     event_trt, event_con, event_trt_long, event_con_long,
                     compare.vec_trt, compare.vec_con){
  Kijv_1 = I_K_ijv/IPCW_t
  temp_sum_1 = NULL
  for(j in 1:num2){
    temp_sum_1 = c(temp_sum_1,sum(I_K_ijv[con_id == j]))
  }
  temp_sum_1_long = temp_sum_1[con_id]
  temp_sum_2_long = compare.vec_trt[trt_id]
  Kijv_2 = ((event_trt_long == 0)*(g1Xc_K >= cen_time_trt_long)*temp_sum_1_long)/
    (IPCW_t * temp_sum_2_long)
  Kijv_3 = 0
  for(p in 1:num1){
    Kijv_3 = Kijv_3 +
      ((event_trt[p] == 0)*(cen_time_trt_long >= cen_time_trt[p])*(g1Xc_K >= cen_time_trt[p])*temp_sum_1_long)/
      (IPCW_t * (compare.vec_trt[p]^2))
  }
  temp_sum_3_long = 0
  for(m in 1:num2){
    temp_sum_3_long = temp_sum_3_long +
      (I_K_ijv[seq(m,num1*num2,num2)][trt_id])*
      (g2Xc_K[m] >= cen_time_con_long)/IPCW_t
  }
  temp_sum_4_long = compare.vec_con[con_id]
  Kijv_4 = ((event_con_long == 0)*temp_sum_3_long)/temp_sum_4_long

  Kijv_5 = 0
  for(q in 1:num2){
    Kijv_5 = Kijv_5 +
      ((event_con[q] == 0)*(cen_time_con_long >= cen_time_con[q])*(temp_sum_3_long[seq(q,num1*num2,num2)][trt_id]))/
      (compare.vec_con[q]^2)
  }

  Kijv = Kijv_1 + Kijv_2 - Kijv_3 + Kijv_4 - Kijv_5
  return(Kijv)
}
