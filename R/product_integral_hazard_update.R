##Some helpful functions
product_integral_func_update<-function(start,end,first,second,third,
                                        basis_list_N3_select_update,coef_N3_update,
                                             A1_time_grid,basis_list_A1_select,coef_A1,
                                             A2_time_grid,basis_list_A2_select,coef_A2,
                                             N1_time_grid,basis_list_N1_select,coef_N1,
                                             N2_time_grid,basis_list_N2_select,coef_N2,
                                       s1_list,s2_list,
                                       tt1,tt2,
                                       delt1,delt2,
                                       grid_scale_clever_covariate,
                                       type,
                                       censoring){
  prob<-1
  if(start<end){
    if(type=="N1"){
      grid_1_N1<-sort(c(N1_time_grid[N1_time_grid<end & N1_time_grid>start],start,end))
      for(i in 1:(length(grid_1_N1)-1)){
        temp<-data.frame(t=(grid_1_N1[i]+grid_1_N1[i+1])/2,tt2_less_t=first,tt2_less_t_delt2=second,tt2_less_t_tt2=third)
        temp2<-hal9001::make_design_matrix(as.matrix(temp),
                                           basis_list_N3_select_update[-c(length(basis_list_N3_select_update):(length(basis_list_N3_select_update)-length(s1_list)+1))])

        clever_covariate_vec<-NULL
        for(j in 1:length(s1_list)){
          clever_covariate<-clever_covaraite_hazards_func(s1=s1_list[[j]],s2=s2_list[[j]],u=(grid_1_N1[i]+grid_1_N1[i+1])/2,tt1=tt1,tt2=tt2,delt1=delt1,delt2=delt2,
                                                          A1_time_grid,basis_list_A1_select,coef_A1,
                                                          A2_time_grid,basis_list_A2_select,coef_A2,
                                                          N1_time_grid,basis_list_N1_select,coef_N1,
                                                          N2_time_grid,basis_list_N2_select,coef_N2,
                                                          censoring=censoring,
                                                          grid_scale_clever_covariate)[1]
          clever_covariate_vec<-c(clever_covariate_vec,clever_covariate)
        }
        temp2<-c(1,as.vector(temp2),clever_covariate_vec)
        log_hazard_midpoint_N1<-sum(coef_N3_update*temp2)
        prob_N1_interval<-exp(-(grid_1_N1[i+1]-grid_1_N1[i])*exp(log_hazard_midpoint_N1))
        prob<-prob*prob_N1_interval
      }
    }else{
      grid_1_N2<-sort(c(N2_time_grid[N2_time_grid<end & N2_time_grid>start],start,end))
      for(i in 1:(length(grid_1_N2)-1)){
        temp<-data.frame(t=(grid_1_N2[i]+grid_1_N2[i+1])/2,tt1_less_t=first,tt1_less_t_delt1=second,tt1_less_t_tt1=third)
        temp2<-hal9001::make_design_matrix(as.matrix(temp),
                                           basis_list_N3_select_update[-c(length(basis_list_N3_select_update):(length(basis_list_N3_select_update)-length(s1_list)+1))])
        clever_covariate_vec<-NULL
        for(j in 1:length(s1_list)){
          clever_covariate<-clever_covaraite_hazards_func(s1=s1_list[[j]],s2=s2_list[[j]],u=(grid_1_N2[i]+grid_1_N2[i+1])/2,tt1=tt1,tt2=tt2,delt1=delt1,delt2=delt2,
                                                          A1_time_grid,basis_list_A1_select,coef_A1,
                                                          A2_time_grid,basis_list_A2_select,coef_A2,
                                                          N1_time_grid,basis_list_N1_select,coef_N1,
                                                          N2_time_grid,basis_list_N2_select,coef_N2,
                                                          censoring=censoring,
                                                          grid_scale_clever_covariate)[2]
          clever_covariate_vec<-c(clever_covariate_vec,clever_covariate)
        }
        temp2<-c(1,as.vector(temp2),clever_covariate_vec)
        log_hazard_midpoint_N2<-sum(coef_N3_update*temp2)
        prob_N2_interval<-exp(-(grid_1_N2[i+1]-grid_1_N2[i])*exp(log_hazard_midpoint_N2))
        prob<-prob*prob_N2_interval
      }
    }
  }
  return(prob)
}


hazard_func_update<-function(t,first,second,third,basis_list_N3_select_update,coef_N3_update,
                             A1_time_grid,basis_list_A1_select,coef_A1,
                             A2_time_grid,basis_list_A2_select,coef_A2,
                             N1_time_grid,basis_list_N1_select,coef_N1,
                             N2_time_grid,basis_list_N2_select,coef_N2,
                             s1_list,s2_list,
                             tt1,tt2,
                             delt1,delt2,
                             grid_scale_clever_covariate,
                             type,
                             censoring){
  if(type=="N1"){
    temp<-data.frame(t=t,tt2_less_t=first,tt2_less_t_delt2=second,tt2_less_t_tt2=third)
  }else{
    temp<-data.frame(t=t,tt1_less_t=first,tt1_less_t_delt1=second,tt1_less_t_tt1=third)
  }

  temp2<-hal9001::make_design_matrix(as.matrix(temp),
                                     basis_list_N3_select_update[-c(length(basis_list_N3_select_update):(length(basis_list_N3_select_update)-length(s1_list)+1))])

  clever_covariate_vec<-NULL
  for(i in 1:length(s1_list)){
    clever_covariates<-clever_covaraite_hazards_func(s1_list[[i]],s2=s2_list[[i]],u=t,tt1=tt1,tt2=tt2,delt1=delt1,delt2=delt2,
                                                    A1_time_grid,basis_list_A1_select,coef_A1,
                                                    A2_time_grid,basis_list_A2_select,coef_A2,
                                                    N1_time_grid,basis_list_N1_select,coef_N1,
                                                    N2_time_grid,basis_list_N2_select,coef_N2,
                                                    censoring=censoring,
                                                    grid_scale_clever_covariate)
    if(type=="N1"){
      clever_covariate<-clever_covariates[1]
    }else{
      clever_covariate<-clever_covariates[2]
    }

    clever_covariate_vec<-c(clever_covariate_vec,clever_covariate)
  }

  temp2<-c(1,as.vector(temp2),clever_covariate_vec)
  log_hazard_N3<-sum(coef_N3_update*temp2)
  hazard_N3<-exp(log_hazard_N3)
  return(hazard_N3)
}




