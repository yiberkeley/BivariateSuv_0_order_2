#Once get the updated conditional intensity estimated, we want to get the updated canonical gradient and the updated estimated survival curve.

Estimation_update<-function(
                      coef_N1_initial_update,basis_list_N1_select_update,
                      coef_N2_initial_update,basis_list_N2_select_update,
                           coef_N1_initial,basis_list_N1_select,
                           coef_N2_initial,basis_list_N2_select,
                               coef_A1_initial,basis_list_A1_select,
                               coef_A2_initial,basis_list_A2_select,
                               s1,s2,
                               data_obs,monitoring_time,
                               grid_scale_clever_covariate
                                #,h_info_old
                          ){


  N1_l1_norm_update<-sum(abs(coef_N1_initial_update[-1]))
  N2_l1_norm_update<-sum(abs(coef_N2_initial_update[-1]))

  N1_num_coef_update<-length(basis_list_N1_select_update)
  N2_num_coef_update<-length(basis_list_N2_select_update)



  #Extact the time points where the fitted conditional intensity is possible to change
  #N1 jumping process
  col_index_temp<-unlist(lapply(basis_list_N1_select, `[[`, 1))
  val_temp<-unlist(lapply(basis_list_N1_select, `[[`, 2))
  N1_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
  #N2 jumping process
  col_index_temp<-unlist(lapply(basis_list_N2_select, `[[`, 1))
  val_temp<-unlist(lapply(basis_list_N2_select, `[[`, 2))
  N2_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
  #A1 jumping process
  col_index_temp<-unlist(lapply(basis_list_A1_select, `[[`, 1))
  val_temp<-unlist(lapply(basis_list_A1_select, `[[`, 2))
  A1_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
  #A2 jumping process
  col_index_temp<-unlist(lapply(basis_list_A2_select, `[[`, 1))
  val_temp<-unlist(lapply(basis_list_A2_select, `[[`, 2))
  A2_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))



  ##indirect approach
  #delt1=0,delt2=0
  density_estimtaed_update<-function(t1,t2,
                              basis_list_N1_select_update,coef_N1_initial_update,
                              basis_list_N2_select_update,coef_N2_initial_update,
                              A1_time_grid,basis_list_A1_select,coef_A1,
                              A2_time_grid,basis_list_A2_select,coef_A2,
                              N1_time_grid,basis_list_N1_select,coef_N1,
                              N2_time_grid,basis_list_N2_select,coef_N2,
                              s1,s2){
    prob_overall<-1

    if(t1<=t2){

      #Before t1
      ##Create a fine grid from 0,t1, including t1
      if(t1!=0){

        prob_overall<-prob_overall*product_integral_func_update(start=0,end=t1,0,0,0,
                                                                basis_list_N1_select_update,coef_N1_initial_update,
                                                                   A1_time_grid,basis_list_A1_select,coef_A1,
                                                                   A2_time_grid,basis_list_A2_select,coef_A2,
                                                                   N1_time_grid,basis_list_N1_select,coef_N1,
                                                                   N2_time_grid,basis_list_N2_select,coef_N2,
                                                                     s1,s2,
                                                                     t1,t2,
                                                                      1,1,
                                                                    type="N1")

        prob_overall<-prob_overall*product_integral_func_update(start=0,end=t1,0,0,0,
                                                                 basis_list_N2_select_update,coef_N2_initial_update,
                                                                 A1_time_grid,basis_list_A1_select,coef_A1,
                                                                 A2_time_grid,basis_list_A2_select,coef_A2,
                                                                 N1_time_grid,basis_list_N1_select,coef_N1,
                                                                 N2_time_grid,basis_list_N2_select,coef_N2,
                                                                 s1,s2,
                                                                 t1,t2,
                                                                  1,1,
                                                                 type="N2")
      }
      #At t1
      prob_overall<-prob_overall*hazard_func_update(t1,0,0,0,basis_list_N1_select_update,coef_N1_initial_update,
                                   A1_time_grid,basis_list_A1_select,coef_A1,
                                   A2_time_grid,basis_list_A2_select,coef_A2,
                                   N1_time_grid,basis_list_N1_select,coef_N1,
                                   N2_time_grid,basis_list_N2_select,coef_N2,
                                   s1,s2,
                                   t1,t2,
                                   1,1,
                                   type="N1")

      #Between t1 and t2
      ##Create a fine grid from t1 to t2
      if(t1<t2){
        prob_overall<-prob_overall*product_integral_func_update(start=t1,end=t2,1,1,t1,
                                                                 basis_list_N2_select_update,coef_N2_initial_update,
                                                                 A1_time_grid,basis_list_A1_select,coef_A1,
                                                                 A2_time_grid,basis_list_A2_select,coef_A2,
                                                                 N1_time_grid,basis_list_N1_select,coef_N1,
                                                                 N2_time_grid,basis_list_N2_select,coef_N2,
                                                                 s1,s2,
                                                                 t1,t2,
                                                                1,1,
                                                                 type="N2")
      }

      #At t2
      prob_overall<-prob_overall*hazard_func_update(t2,1,1,t1,basis_list_N2_select_update,coef_N2_initial_update,
                                                     A1_time_grid,basis_list_A1_select,coef_A1,
                                                     A2_time_grid,basis_list_A2_select,coef_A2,
                                                     N1_time_grid,basis_list_N1_select,coef_N1,
                                                     N2_time_grid,basis_list_N2_select,coef_N2,
                                                     s1,s2,
                                                     t1,t2,
                                                      1,1,
                                                     type="N2")

    }else{
      #Before t2
      ##Create a fine grid from 0,t2, including t2
      if(t2!=0){
        prob_overall<-prob_overall*product_integral_func_update(start=0,end=t2,0,0,0,
                                                                basis_list_N1_select_update,coef_N1_initial_update,
                                                                A1_time_grid,basis_list_A1_select,coef_A1,
                                                                A2_time_grid,basis_list_A2_select,coef_A2,
                                                                N1_time_grid,basis_list_N1_select,coef_N1,
                                                                N2_time_grid,basis_list_N2_select,coef_N2,
                                                                s1,s2,
                                                                t1,t2,
                                                                1,1,
                                                                type="N1")

        prob_overall<-prob_overall*product_integral_func_update(start=0,end=t2,0,0,0,
                                                                basis_list_N2_select_update,coef_N2_initial_update,
                                                                A1_time_grid,basis_list_A1_select,coef_A1,
                                                                A2_time_grid,basis_list_A2_select,coef_A2,
                                                                N1_time_grid,basis_list_N1_select,coef_N1,
                                                                N2_time_grid,basis_list_N2_select,coef_N2,
                                                                s1,s2,
                                                                t1,t2,
                                                                1,1,
                                                                type="N2")
      }
      #At t2
      prob_overall<-prob_overall*hazard_func_update(t2,0,0,0,basis_list_N2_select_update,coef_N2_initial_update,
                                                     A1_time_grid,basis_list_A1_select,coef_A1,
                                                     A2_time_grid,basis_list_A2_select,coef_A2,
                                                     N1_time_grid,basis_list_N1_select,coef_N1,
                                                     N2_time_grid,basis_list_N2_select,coef_N2,
                                                     s1,s2,
                                                     t1,t2,
                                                      1,1,
                                                     type="N2")

      #Between t2 and t1
      ##Create a fine grid from t2 to t1
      prob_overall<-prob_overall*product_integral_func_update(start=t2,end=t1,1,1,t2,
                                                                 basis_list_N1_select_update,coef_N1_initial_update,
                                                                 A1_time_grid,basis_list_A1_select,coef_A1,
                                                                 A2_time_grid,basis_list_A2_select,coef_A2,
                                                                 N1_time_grid,basis_list_N1_select,coef_N1,
                                                                 N2_time_grid,basis_list_N2_select,coef_N2,
                                                                 s1,s2,
                                                                 t1,t2,
                                                                 1,1,
                                                                 type="N1")

      #At t1
      prob_overall<-prob_overall*hazard_func_update(t1,1,1,t2,basis_list_N1_select_update,coef_N1_initial_update,
                                                       A1_time_grid,basis_list_A1_select,coef_A1,
                                                       A2_time_grid,basis_list_A2_select,coef_A2,
                                                       N1_time_grid,basis_list_N1_select,coef_N1,
                                                       N2_time_grid,basis_list_N2_select,coef_N2,
                                                       s1,s2,
                                                       t1,t2,
                                                        1,1,
                                                       type="N1")

    }
    return(prob_overall)
  }


  #Get the densities
  #This is the grid for the density grid
  grid_scale_density<-0.01

  t1_density_grid<- seq(0,1,grid_scale_density)
  t2_density_grid <- seq(0,1,grid_scale_density)

  density_data<-expand.grid(t1=t1_density_grid,t2=t2_density_grid)


  system.time(

    density_grid<-future_apply(density_data,1,function(x)density_estimtaed_update(t1=as.numeric(x[1]),t2=as.numeric(x[2]),
                                                                                  basis_list_N1_select_update,coef_N1_initial_update,
                                                                                  basis_list_N2_select_update,coef_N2_initial_update,
                                                                                  A1_time_grid,basis_list_A1_select,coef_A1_initial,
                                                                                  A2_time_grid,basis_list_A2_select,coef_A2_initial,
                                                                                  N1_time_grid,basis_list_N1_select,coef_N1_initial,
                                                                                  N2_time_grid,basis_list_N2_select,coef_N2_initial,
                                                                                  s1,s2))
  )

  density_data_complete<-cbind(density_data,density_grid)

  #Visualization
  #plotly::plot_ly(x=density_data_complete$t1,y=density_data_complete$t2,z=density_data_complete$density_grid)
  #With library(ggfx)
  # ggplot(density_data_complete, aes(x=t1, y=t2))+with_blur(
  #   geom_raster(aes(fill = density_grid), interpolate = FALSE),
  #   sigma = 10
  # ) +scale_fill_viridis_c(option = "inferno")

  #Survival function through the density
  #left riemman sum

  S_estimate_2<-function(t_1,t_2){
    temp_data<-density_data_complete%>%filter(t1>t_1,t2>t_2)
    return(sum(temp_data$density_grid)*grid_scale_density*grid_scale_density)
  }
  #right riemman sum
  S_estimate_3<-function(t_1,t_2){
    temp_data<-density_data_complete%>%filter(t1>=t_1,t2>=t_2,t1<1,t2<1)
    return(sum(temp_data$density_grid)*grid_scale_density*grid_scale_density)
  }

  #Plot the fitted survival curve -Indirect Approach
  t1_plot<- seq(0,1,0.1)
  t2_plot <- seq(0,1,0.1)

  plot_data<-expand.grid(t1=t1_plot,t2=t2_plot)


  plan(multisession,workers = nCores)
  survival_curve_2<-future_apply(plot_data,1,function(x)S_estimate_2(x[1],x[2]))


  surve_indirect<-expand.grid(t1=t1_plot,t2=t2_plot)
  surve_indirect<-surve_indirect%>%mutate(survival_curve=survival_curve_2)
  #Visualization
  # plotly::plot_ly(x=surve_indirect$t1,y=surve_indirect$t2,z=surve_indirect$survival_curve)
  # ggplot2::ggplot(surve_indirect, aes(x=t1, y=t2))+with_blur(
  #   geom_raster(aes(fill = survival_curve), interpolate = FALSE),
  #   sigma = 10
  # ) +scale_fill_viridis_c(option = "inferno")


  ##Now, get the clever covariate at s1,s2 of interest and also the D* info


  h_info_update<-H_Dstar_info_func_update(s1,s2,data_obs,monitoring_time,
                                     A1_time_grid,basis_list_A1_select,coef_A1_initial,
                                     A2_time_grid,basis_list_A2_select,coef_A2_initial,
                                     N1_time_grid,basis_list_N1_select,coef_N1_initial,
                                     N2_time_grid,basis_list_N2_select,coef_N2_initial,
                                     grid_scale_clever_covariate,
                                     basis_list_N1_select_update,coef_N1_initial_update,
                                     basis_list_N2_select_update,coef_N2_initial_update)
                                     #,h_info_old)




  result_inner<-list()
  result_inner[[1]]<-surve_indirect
  result_inner[[2]]<-N1_l1_norm_update
  result_inner[[3]]<-N2_l1_norm_update
  result_inner[[4]]<-N1_num_coef_update
  result_inner[[5]]<-N2_num_coef_update
  result_inner[[6]]<-h_info_update

  # result_inner[[6]]<-cv_selected_cv_risk_N1
  # result_inner[[7]]<-cv_selected_cv_risk_N2
  # result_inner[[8]]<-cv_selected_cv_number_coef_N1
  # result_inner[[9]]<-cv_selected_cv_number_coef_N2
  # result_inner[[10]]<-cv_selected_fit_risk_N1
  # result_inner[[11]]<-cv_selected_fit_risk_N2



  #print("one survival fit done in undersmoothing")


  return(result_inner)
}
