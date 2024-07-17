#Once get the conditional intensity estimated, we want to get the canonical gradient and the estimated survival curve.

Estimation<-function(coef_N1_initial,basis_list_N1_select,
                     coef_N2_initial,basis_list_N2_select,
                     coef_A1_initial,basis_list_A1_select,
                     coef_A2_initial,basis_list_A2_select,
                     censoring,
                     s1_list,
                     s2_list,
                     data_obs,monitoring_time,
                     grid_scale_clever_covariate,
                     targeting){


        N1_l1_norm<-sum(abs(coef_N1_initial[-1]))
        N2_l1_norm<-sum(abs(coef_N2_initial[-1]))

        N1_num_coef<-length(basis_list_N1_select)
        N2_num_coef<-length(basis_list_N2_select)



        #Extact the time points where the fitted conditional intensity is possible to change
        #N1 jumping process
        col_index_temp<-unlist(lapply(basis_list_N1_select, `[[`, 1))
        val_temp<-unlist(lapply(basis_list_N1_select, `[[`, 2))
        N1_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
        #N2 jumping process
        col_index_temp<-unlist(lapply(basis_list_N2_select, `[[`, 1))
        val_temp<-unlist(lapply(basis_list_N2_select, `[[`, 2))
        N2_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
        if(censoring){
          #A1 jumping process
          col_index_temp<-unlist(lapply(basis_list_A1_select, `[[`, 1))
          val_temp<-unlist(lapply(basis_list_A1_select, `[[`, 2))
          A1_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
          #A2 jumping process
          col_index_temp<-unlist(lapply(basis_list_A2_select, `[[`, 1))
          val_temp<-unlist(lapply(basis_list_A2_select, `[[`, 2))
          A2_time_grid<-sort(unique(c(0,val_temp[col_index_temp==1|col_index_temp==4])))
        }else{
          A1_time_grid<-NULL
          A2_time_grid<-NULL
        }



        ##indirect approach
        #delt1=0,delt2=0
        density_estimtaed<-function(t1,t2,coef_N1,coef_N2,basis_list_N1_select,basis_list_N2_select,N1_time_grid,N2_time_grid){
          prob_overall<-1
          if(t1<=t2){

            #Before t1
            ##Create a fine grid from 0,t1, including t1
            if(t1!=0){
              prob_overall<-prob_overall*product_integral_func(start=0,end=t1,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
              prob_overall<-prob_overall*product_integral_func(start=0,end=t1,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
            }
            #At t1
            prob_overall<-prob_overall*hazard_func(t1,0,0,0,basis_list_N1_select,coef_N1,type="N1")

            #Between t1 and t2
            ##Create a fine grid from t1 to t2
            if(t1<t2){
              prob_overall<-prob_overall*product_integral_func(start=t1,end=t2,1,1,t1,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
              #At t2
              prob_overall<-prob_overall*hazard_func(t2,1,1,t1,basis_list_N2_select,coef_N2,type="N2")

            }else{
              #################
              prob_overall<-prob_overall*hazard_func(t2,0,0,0,basis_list_N2_select,coef_N2,type="N2")
              #################
            }



          }else{
            #Before t2
            ##Create a fine grid from 0,t2, including t2
            if(t2!=0){
              prob_overall<-prob_overall*product_integral_func(start=0,end=t2,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")
              prob_overall<-prob_overall*product_integral_func(start=0,end=t2,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type="N2")
            }
            #At t2
            prob_overall<-prob_overall*hazard_func(t2,0,0,0,basis_list_N2_select,coef_N2,type="N2")

            #Between t2 and t1
            ##Create a fine grid from t2 to t1
            prob_overall<-prob_overall*product_integral_func(start=t2,end=t1,1,1,t2,N1_time_grid,basis_list_N1_select,coef_N1,type="N1")

            #At t1
            prob_overall<-prob_overall*hazard_func(t1,1,1,t2,basis_list_N1_select,coef_N1,type="N1")

          }

          #check point
          if(t1==0.1 & t2==0.1){
            print("check point: t1=0.1 t2=0.1 initial density fit finished")
          }
          if(t1==0.5 & t2==0.5){
            print("check point: t1=0.5 t2=0.5 initial density fit finished")
          }
          if(t1==0.9 & t2==0.9){
            print("check point: t1=0.9 t2=0.9 initial density fit finished")
          }

          return(prob_overall)
        }


        #Get the densities
        #This is the grid for the density grid
        grid_scale_density<-0.01 #used to be 0.01

        t1_density_grid<- seq(0,1,grid_scale_density)[-1]
        t2_density_grid <- seq(0,1,grid_scale_density)[-1]

        density_data<-expand.grid(t1=t1_density_grid,t2=t2_density_grid)

        #nCores<-detectCores()-2
        #plan(multisession,workers = nCores)

        print("check point: prepare to fit initial densities")

        plan(multicore)

        system.time(

          density_grid<-future.apply::future_apply(density_data,1,function(x)density_estimtaed(t1=as.numeric(x[1]),t2=as.numeric(x[2]),
                                                                                       coef_N1=coef_N1_initial,
                                                                                       coef_N2=coef_N2_initial,
                                                                                       basis_list_N1_select,
                                                                                       basis_list_N2_select,
                                                                                       N1_time_grid,
                                                                                       N2_time_grid))
        )

        density_data_complete<-cbind(density_data,density_grid)

        print("check point: initial density fit finished")

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
          temp_data<-density_data_complete%>%filter(t1>t_1,t2>t_2)%>%mutate(length1=grid_scale_density,
                                                                            length2=grid_scale_density)
          return(sum(temp_data$density_grid*temp_data$length1*temp_data$length2))
        }

              #Also report the CDF estimate
              CDF_estimate_2<-function(t_1,t_2){
                temp_data<-density_data_complete%>%filter(0<t1,t1<=t_1,0<t2,t2<=t_2)%>%mutate(length1=grid_scale_density,
                                                                                  length2=grid_scale_density)
                return(sum(temp_data$density_grid*temp_data$length1*temp_data$length2))
              }

        #right riemman sum
        S_estimate_3<-function(t_1,t_2){
          temp_data<-density_data_complete%>%filter(t1>=t_1,t2>=t_2,t1<1,t2<1)
          return(sum(temp_data$density_grid)*grid_scale_density*grid_scale_density)
        }

        #Plot the fitted survival curve -Indirect Approach
        t1_plot<- seq(0,1,0.01)
        t2_plot <- seq(0,1,0.01)

        plot_data<-expand.grid(t1=t1_plot,t2=t2_plot)

        #nCores<-detectCores()-2
        #plan(multisession,workers = nCores)

        plan(multicore)
        survival_curve_2<-future.apply::future_apply(plot_data,1,function(x)S_estimate_2(x[1],x[2]))

        surve_indirect<-expand.grid(t1=t1_plot,t2=t2_plot)
        surve_indirect<-surve_indirect%>%mutate(survival_curve=survival_curve_2)

            plan(multicore)
            cdf_curve_2<-future.apply::future_apply(plot_data,1,function(x)CDF_estimate_2(x[1],x[2]))

            cdf_indirect<-expand.grid(t1=t1_plot,t2=t2_plot)
            cdf_indirect<-cdf_indirect%>%mutate(cdf_curve=cdf_curve_2)


        print("check point: initial survival fit finished")

        #Visualization
        #plotly::plot_ly(x=surve_indirect$t1,y=surve_indirect$t2,z=surve_indirect$survival_curve)
        # ggplot(surve_indirect, aes(x=t1, y=t2))+with_blur(
        #   geom_raster(aes(fill = survival_curve), interpolate = FALSE),
        #   sigma = 10
        # ) +scale_fill_viridis_c(option = "inferno")

        ##################################################
        #########################
        if(targeting){

          ##Now, get the clever covariate at s1,s2 of interest and also the D* info
          `%+%` <- function(a, b) paste0(a, b)

          h_info_list<-list()

          for(i in 1:length(s1_list)){
            print("check piont: prepare to get clever covariates for point " %+% i)
            h_info<-H_Dstar_info_func(s1_list[[i]],s2_list[[i]],
                                      data_obs,monitoring_time,
                                      A1_time_grid,basis_list_A1_select,coef_A1_initial,
                                      A2_time_grid,basis_list_A2_select,coef_A2_initial,
                                      N1_time_grid,basis_list_N1_select,coef_N1_initial,
                                      N2_time_grid,basis_list_N2_select,coef_N2_initial,
                                      censoring=censoring,
                                      grid_scale_clever_covariate)

            h_info_list[[i]]<-h_info
            print("check piont: clever covariates done for point " %+% i)
          }



          result_inner<-list()
          result_inner[[1]]<-surve_indirect
          result_inner[[2]]<-N1_l1_norm
          result_inner[[3]]<-N2_l1_norm
          result_inner[[4]]<-N1_num_coef
          result_inner[[5]]<-N2_num_coef
          result_inner[[6]]<-h_info_list
          result_inner[[7]]<-cdf_indirect
          result_inner[[8]]<-density_data_complete


        }else{
          result_inner<-list()
          result_inner[[1]]<-surve_indirect
          result_inner[[2]]<-N1_l1_norm
          result_inner[[3]]<-N2_l1_norm
          result_inner[[4]]<-N1_num_coef
          result_inner[[5]]<-N2_num_coef
          result_inner[[6]]<-(-1)
          result_inner[[7]]<-cdf_indirect
          result_inner[[8]]<-density_data_complete
        }
        #########################
        ##################################################

        # result_inner<-list()
        # result_inner[[1]]<-surve_indirect
        # result_inner[[2]]<-N1_l1_norm
        # result_inner[[3]]<-N2_l1_norm
        # result_inner[[4]]<-N1_num_coef
        # result_inner[[5]]<-N2_num_coef
        # result_inner[[6]]<-(-1)
        # result_inner[[7]]<-cdf_indirect
        # result_inner[[8]]<-density_data_complete


      return(result_inner)
}


