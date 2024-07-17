
#Update using a densor sampling of monitoring time when do the update.


Targeting<-function(data_obs,foldid,perc_sample,monitoring_time,
                    basis_list_N1_select,basis_list_N2_select,
                    coef_N1_initial,coef_N2_initial,
                    penalty_N1_select,penalty_N2_select,
                    h_info_list,
                    undersmooth_type="cv_all",
                    cv_choice="lambda",
                    nfolds,
                    masking_option,
                    weight_option,
                    cvxr,
                    l1_norm_initial_N1,
                    l1_norm_initial_N2,
                    uniform_ratio,
                    relaxed,
                    relaxed_initial=F,
                    offset_initial=NULL,
                    offset_relax=NULL){

            poisson_cv_fit_func_all<-function(cvxr,poisson_data_N3,basis_list_N3,nfolds,penalty.factor_N3,undersmooth=F,l1_norm_start=F,uniform_ratio){
              if(!cvxr){
                result<-poisson_cv_fit_func(poisson_data_N3,basis_list_N3,nfolds,penalty.factor_N3,undersmooth,l1_norm_start,uniform_ratio)
              }else{
                result<-poisson_cv_fit_func_cvxr(poisson_data_N3,basis_list_N3,nfolds,penalty.factor_N3,undersmooth,l1_norm_start,uniform_ratio)
              }
              return(result)
            }

          ##Create repeated data
          repeated_data_update_collection<-repeated_data_func(data_obs,foldid,min(perc_sample*9,1),monitoring_time,censor_data_include=F,masking_option,weight_option)
          repeated_data_N1_2_update<-repeated_data_update_collection[[1]]
          repeated_data_N2_2_update<-repeated_data_update_collection[[2]]
          #The only possible error is coming from the fact that it t2 is btw two monitoring time sampled and thus there is some approximation
          #in that interval for the hazard calculation in that interval, but we should be fine.



          ##New poisson data
          #Now,augment the Repeated data with basis selected and the clever covariate calculated
          #Poisson data N1
          basis_expansion<-make_design_matrix(as.matrix(repeated_data_N1_2_update%>%select(t,tt2_less_t,tt2_less_t_delt2,tt2_less_t_tt2)),basis_list_N1_select)
          basic_info<-repeated_data_N1_2_update
          poisson_data_N1_update<-cbind(basic_info,as.matrix(basis_expansion))
          #Poisson data N2
          basis_expansion<-make_design_matrix(as.matrix(repeated_data_N2_2_update%>%select(t,tt1_less_t,tt1_less_t_delt1,tt1_less_t_tt1)),basis_list_N2_select)
          basic_info<-repeated_data_N2_2_update
          poisson_data_N2_update<-cbind(basic_info,as.matrix(basis_expansion))

          if(!relaxed_initial){
              #Clever covariates info
              clever_covariate_table_list<-list()

              for(i in 1:length(h_info_list)){

                h_info<-h_info_list[[i]]

                obs_id_vec<-lapply(h_info, `[[`, 1)
                obs_id_vec<-as.vector(unlist(obs_id_vec))

                time_point_vec<-lapply(h_info, `[[`, 2)
                time_point_vec<-as.vector(unlist(time_point_vec))

                H1_vec<-lapply(h_info, `[[`, 3)
                H1_vec<-as.vector(unlist(H1_vec))

                H2_vec<-lapply(h_info, `[[`, 4)
                H2_vec<-as.vector(unlist(H2_vec))

                clever_covariate_table<-data.frame(obs=obs_id_vec,t=time_point_vec,H1=H1_vec,H2=H2_vec)

                clever_covariate_table_list[[i]]<-clever_covariate_table
              }


              #add clever covariate info to poisson data
              #Poisson data N1

              poisson_data_N1_update_clever<-poisson_data_N1_update
              for(i in 1:length(h_info_list)){
                poisson_data_N1_update_clever <- merge(x=poisson_data_N1_update_clever,y=clever_covariate_table_list[[i]],
                                                       by=c("obs","t"),
                                                       all.x=TRUE)
                poisson_data_N1_update_clever[is.na(poisson_data_N1_update_clever)]<-0
                poisson_data_N1_update_clever<-poisson_data_N1_update_clever%>%select(-H2)
                poisson_data_N1_update_clever<-poisson_data_N1_update_clever%>%rename(!! paste0("H1.",i):=H1)
              }

              #Poisson data N2
              poisson_data_N2_update_clever<-poisson_data_N2_update
              for(i in 1:length(h_info_list)){
                  poisson_data_N2_update_clever <- merge(x=poisson_data_N2_update_clever,y=clever_covariate_table_list[[i]],
                                                         by=c("obs","t"),
                                                         all.x=TRUE)
                  poisson_data_N2_update_clever[is.na(poisson_data_N2_update_clever)]<-0
                  poisson_data_N2_update_clever<-poisson_data_N2_update_clever%>%select(-H1)
                  poisson_data_N2_update_clever<-poisson_data_N2_update_clever%>%rename(!! paste0("H2.",i):=H2)
              }



              ##Set all penalty to be 1 except the clever covariate to be 0
              #For numerical issues, we penalize the clever covairates as 0.01
              penalty.factor_N1_update<-c(penalty_N1_select, rep(0,length(h_info_list))) #c(0.005,0.03,0.03,0.06))
              penalty.factor_N2_update<-c(penalty_N2_select, rep(0,length(h_info_list))) #c(0.005,0.03,0.03,0.06))


              H1_names<-names(poisson_data_N1_update_clever)[grepl(pattern = "^H1", names(poisson_data_N1_update_clever))]
              H2_names<-names(poisson_data_N2_update_clever)[grepl(pattern = "^H2", names(poisson_data_N2_update_clever))]

          }else{
              poisson_data_N1_update_clever<-poisson_data_N1_update
              poisson_data_N2_update_clever<-poisson_data_N2_update

              penalty.factor_N1_update<-penalty_N1_select
              penalty.factor_N2_update<-penalty_N2_select

              H1_names<-NULL
              H2_names<-NULL
          }
          print("check point: Targeting step is prepared")

          ##Targeting
          if(undersmooth_type=="cv_all"){
            if(cv_choice=="lambda"){
              #N1 jumping process
              if(!relaxed & !offset_initial){
                  #CV version of targeting
                  result_temp<-poisson_cv_fit_func_all(cvxr,poisson_data_N1_update_clever,c(basis_list_N1_select,H1_names),nfolds,penalty.factor_N1_update,undersmooth=F,l1_norm_start=l1_norm_initial_N1,uniform_ratio)
                  cv_selected_cv_risk_N1_update <- result_temp[[1]]
                  cv_selected_cv_number_coef_N1_update <- result_temp[[2]]
                  cv_selected_fit_risk_N1_update <- result_temp[[3]]
                  coef_N1_initial_list_update <- result_temp[[4]]
                  basis_list_N1_select_list_update <- result_temp[[5]]
              }else if(offset_initial){

                  print("****enter the offset with offset relax value as " %+% offset_relax)
                  result_temp<-poisson_cv_fit_func_cvxr(poisson_data_N3=poisson_data_N1_update_clever,
                                                        basis_list_N3=c(basis_list_N1_select,H1_names),
                                                        nfolds=nfolds,
                                                        penalty.factor_N3=penalty.factor_N1_update,
                                                        undersmooth=F,
                                                        l1_norm_start=F,
                                                        uniform_ratio=uniform_ratio,
                                                        offset_initial=offset_initial,
                                                        length_clever=length(h_info_list),
                                                        coef_initial=coef_N1_initial,
                                                        offset_relax=offset_relax)

                  cv_selected_cv_risk_N1_update <- result_temp[[1]]
                  cv_selected_cv_number_coef_N1_update <- result_temp[[2]]
                  cv_selected_fit_risk_N1_update <- result_temp[[3]]
                  coef_N1_initial_list_update <- result_temp[[4]]
                  basis_list_N1_select_list_update <- result_temp[[5]]

              }else{
                  #Relax version of targeting
                  offset<-log(poisson_data_N1_update_clever$timeIntervalLength)
                  weights<- poisson_data_N1_update_clever$weights
                  y<-poisson_data_N1_update_clever$jump
                  x_basis_poisson<-as.matrix(poisson_data_N1_update_clever[,-c(1:9)])
                  time_vec<-as.vector(poisson_data_N1_update_clever$t)

                  fit_N3<-cvrx_poisson_l1(l1_norm_seq=l1_norm_initial_N1*20,x_basis_poisson,offset,weights,y,penalty.factor_N1_update,time_vec,uniform_ratio)


                  beta<-fit_N3[[1]][[1]]
                  beta_0<-fit_N3[[1]][[2]]
                  cv_selected_fit_risk_N3<-fit_N3[[1]][[3]]

                  coef_N3_initial_list<-list()
                  basis_list_N3_select_list<-list()
                  penalty_N3_select_list<-list()
                  for(i in 1:length(l1_norm_initial_N1)){
                    coef_N3_initial_list[[i]]<-c(beta_0,beta[abs(beta)>10^-4])
                    basis_list_N3_select_list[[i]]<-c(basis_list_N1_select,H1_names)[abs(beta)>10^-4]
                    penalty_N3_select_list[[i]]<-penalty.factor_N1_update[abs(beta)>10^-4]
                  }


                  cv_selected_cv_risk_N1_update <- 1
                  cv_selected_cv_number_coef_N1_update <- 1
                  cv_selected_fit_risk_N1_update <- cv_selected_fit_risk_N3
                  coef_N1_initial_list_update <- coef_N3_initial_list
                  basis_list_N1_select_list_update <- basis_list_N3_select_list
              }
              print("N1 update fit finished")

              #N2 jumping process
              if(!relaxed & !offset_initial){
                  #CV version of targeting
                  result_temp<-poisson_cv_fit_func_all(cvxr,poisson_data_N2_update_clever,c(basis_list_N2_select,H2_names),nfolds,penalty.factor_N2_update,undersmooth=F,l1_norm_start=l1_norm_initial_N2,uniform_ratio)
                  cv_selected_cv_risk_N2_update <- result_temp[[1]]
                  cv_selected_cv_number_coef_N2_update <- result_temp[[2]]
                  cv_selected_fit_risk_N2_update <- result_temp[[3]]
                  coef_N2_initial_list_update <- result_temp[[4]]
                  basis_list_N2_select_list_update <- result_temp[[5]]
              }else if(offset_initial){
                length_clever<-length(h_info_list)
                coef_initial<-coef_N2_initial
                result_temp<-poisson_cv_fit_func_cvxr(poisson_data_N2_update_clever,
                                                      c(basis_list_N2_select,H2_names),
                                                      nfolds,penalty.factor_N2_update,
                                                      undersmooth=F,
                                                      l1_norm_start=F,
                                                      uniform_ratio,
                                                      offset_initial,
                                                      length_clever,
                                                      coef_initial,
                                                      offset_relax=offset_relax)
                cv_selected_cv_risk_N2_update <- result_temp[[1]]
                cv_selected_cv_number_coef_N2_update <- result_temp[[2]]
                cv_selected_fit_risk_N2_update <- result_temp[[3]]
                coef_N2_initial_list_update <- result_temp[[4]]
                basis_list_N2_select_list_update <- result_temp[[5]]

              }else{
                  offset<-log(poisson_data_N2_update_clever$timeIntervalLength)
                  weights<- poisson_data_N2_update_clever$weights
                  y<-poisson_data_N2_update_clever$jump
                  x_basis_poisson<-as.matrix(poisson_data_N2_update_clever[,-c(1:9)])
                  time_vec<-as.vector(poisson_data_N2_update_clever$t)

                  fit_N3<-cvrx_poisson_l1(l1_norm_seq=l1_norm_initial_N2*20,x_basis_poisson,offset,weights,y,penalty.factor_N2_update,time_vec,uniform_ratio)


                  beta<-fit_N3[[1]][[1]]
                  beta_0<-fit_N3[[1]][[2]]
                  cv_selected_fit_risk_N3<-fit_N3[[1]][[3]]

                  coef_N3_initial_list<-list()
                  basis_list_N3_select_list<-list()
                  penalty_N3_select_list<-list()
                  for(i in 1:length(l1_norm_initial_N2)){
                    coef_N3_initial_list[[i]]<-c(beta_0,beta[abs(beta)>10^-4])
                    basis_list_N3_select_list[[i]]<-c(basis_list_N2_select,H2_names)[abs(beta)>10^-4]
                    penalty_N3_select_list[[i]]<-penalty.factor_N2_update[abs(beta)>10^-4]
                  }


                  cv_selected_cv_risk_N2_update <- 1
                  cv_selected_cv_number_coef_N2_update <- 1
                  cv_selected_fit_risk_N2_update <- cv_selected_fit_risk_N3
                  coef_N2_initial_list_update <- coef_N3_initial_list
                  basis_list_N2_select_list_update <- basis_list_N3_select_list
              }
              print("N2 update fit finished")
            }
          }

        result<-list()
        result[[1]]<-coef_N1_initial_list_update
        result[[2]]<-basis_list_N1_select_list_update
        result[[3]]<-coef_N2_initial_list_update
        result[[4]]<-basis_list_N2_select_list_update

        return(result)
}

