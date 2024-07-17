
##Canonical Gradient Evaluation
#s1,s2 are the point of interest for Canonical Gradient
#tt1,tt2,delt1,delt2 are observation i



H_Dstar_individual_info_func_update<-function(s1,s2,tt1,tt2,delt1,delt2,obs_id,monitoring_time,
                                       A1_time_grid,basis_list_A1_select,coef_A1,
                                       A2_time_grid,basis_list_A2_select,coef_A2,
                                       N1_time_grid,basis_list_N1_select,coef_N1,
                                       N2_time_grid,basis_list_N2_select,coef_N2,
                                       grid_scale_clever_covariate,
                                       basis_list_N1_select_update,coef_N1_initial_update,
                                       basis_list_N2_select_update,coef_N2_initial_update
){
  #we only need to evaluate H from [0,max{tt1,tt2, t1-,t2-}]
  #if our monitoring_time is replaced later by a sub sample of it, we need to ensure the 0 and tt1,tt2 are there
  grid <-sort(unique((c(monitoring_time[monitoring_time<=max(s1,s2)],tt1,tt2,0))))

  new_grid<-grid[grid<=max(s1,s2)]
  system.time(
    H_table<-sapply(new_grid,function(x){

      clever_covaraite_hazards_func_update (s1,s2,x,tt1,tt2,delt1,delt2,
                                     A1_time_grid,basis_list_A1_select,coef_A1,
                                     A2_time_grid,basis_list_A2_select,coef_A2,
                                     N1_time_grid,basis_list_N1_select,coef_N1,
                                     N2_time_grid,basis_list_N2_select,coef_N2,
                                     grid_scale_clever_covariate,
                                     basis_list_N1_select_update,coef_N1_initial_update,
                                     basis_list_N2_select_update,coef_N2_initial_update
      )

    })
  )


  canonical_gradient<-sum((H_table[1,]*H_table[3,])*diff(c(new_grid,max(s1,s2))))+
    sum((H_table[2,]*H_table[4,])*diff(c(new_grid,max(s1,s2))))+
    ifelse(delt1==1 & tt1<=max(s1,s2), H_table[3,][which(new_grid==tt1)] ,0)+
    ifelse(delt2==1 & tt2<=max(s1,s2), H_table[4,][which(new_grid==tt2)] ,0)


  result_temp<-list()
  result_temp[[1]]<-rep(obs_id,length(new_grid))
  result_temp[[2]]<-new_grid
  result_temp[[3]]<-H_table[1,]
  result_temp[[4]]<-H_table[2,]
  result_temp[[5]]<-canonical_gradient
  return(result_temp)
}



H_Dstar_info_func_update<-function(s1,s2,data_obs,monitoring_time,
                            A1_time_grid,basis_list_A1_select,coef_A1,
                            A2_time_grid,basis_list_A2_select,coef_A2,
                            N1_time_grid,basis_list_N1_select,coef_N1,
                            N2_time_grid,basis_list_N2_select,coef_N2,
                            grid_scale_clever_covariate,
                            basis_list_N1_select_update,coef_N1_initial_update,
                            basis_list_N2_select_update,coef_N2_initial_update){



  data_obs_temp<-data_obs
  data_obs_temp$obs_id<-1:dim(data_obs)[1]

  system.time(
    h_info<-future_apply(data_obs_temp,1,function(x){
      H_Dstar_individual_info_func_update(s1,s2,tt1=x[1],tt2=x[3],delt1=x[2],delt2=x[4],
                                   obs_id=x[5],
                                   monitoring_time,
                                   A1_time_grid,basis_list_A1_select,coef_A1,
                                   A2_time_grid,basis_list_A2_select,coef_A2,
                                   N1_time_grid,basis_list_N1_select,coef_N1,
                                   N2_time_grid,basis_list_N2_select,coef_N2,
                                   grid_scale_clever_covariate,
                                   basis_list_N1_select_update,coef_N1_initial_update,
                                   basis_list_N2_select_update,coef_N2_initial_update
      )
    })
  )

  return(h_info)
}




