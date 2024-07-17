direct_canonical_update_eval<-function(data_obs,foldid,monitoring_time,
                   coef_N1_initial_update,basis_list_N1_select_update,
                   coef_N2_initial_update,basis_list_N2_select_update,
                     h_info_old_list,
                   masking_option,
                   weight_option){
  ##Create repeated data
  repeated_data_update_collection<-repeated_data_func(data_obs,foldid,perc_sample=1,monitoring_time,censor_data_include=F,masking_option,weight_option)
  repeated_data_N1_2_update<-repeated_data_update_collection[[1]]
  repeated_data_N2_2_update<-repeated_data_update_collection[[2]]


  #The only possible error is coming from the fact that it t2 is btw two monitoring time sampled and thus there is some approximation
  #in that interval for the hazard calculation in that interval, but we should be fine.



  ##New poisson data
  #Now,augment the Repeated data with basis selected and the clever covariate calculated
  #Poisson data N1
  basis_expansion<-make_design_matrix(as.matrix(repeated_data_N1_2_update%>%select(t,tt2_less_t,tt2_less_t_delt2,tt2_less_t_tt2)),
                                      basis_list_N1_select_update[-c(length(basis_list_N1_select_update):(length(basis_list_N1_select_update)-length(h_info_old_list)+1))])
  basic_info<-repeated_data_N1_2_update
  poisson_data_N1_update<-cbind(basic_info,as.matrix(basis_expansion))
  #Poisson data N2
  basis_expansion<-make_design_matrix(as.matrix(repeated_data_N2_2_update%>%select(t,tt1_less_t,tt1_less_t_delt1,tt1_less_t_tt1)),
                                      basis_list_N2_select_update[-c(length(basis_list_N2_select_update):(length(basis_list_N2_select_update)-length(h_info_old_list)+1))])
  basic_info<-repeated_data_N2_2_update
  poisson_data_N2_update<-cbind(basic_info,as.matrix(basis_expansion))

  ##Get the full poisson clever first
            #Clever covariates info
            clever_covariate_table_list<-list()

            for(i in 1:length(h_info_old_list)){

              h_info<-h_info_old_list[[i]]

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

            poisson_data_N1_update_clever_full<-poisson_data_N1_update
            for(i in 1:length(h_info_old_list)){
              poisson_data_N1_update_clever_full <- merge(x=poisson_data_N1_update_clever_full,y=clever_covariate_table_list[[i]],
                                                     by=c("obs","t"),
                                                     all.x=TRUE)
              poisson_data_N1_update_clever_full[is.na(poisson_data_N1_update_clever_full)]<-0
              poisson_data_N1_update_clever_full<-poisson_data_N1_update_clever_full%>%select(-H2)
              poisson_data_N1_update_clever_full<-poisson_data_N1_update_clever_full%>%rename(!! paste0("H1.",i):=H1)
            }

            #Poisson data N2
            poisson_data_N2_update_clever_full<-poisson_data_N2_update
            for(i in 1:length(h_info_old_list)){
              poisson_data_N2_update_clever_full <- merge(x=poisson_data_N2_update_clever_full,y=clever_covariate_table_list[[i]],
                                                     by=c("obs","t"),
                                                     all.x=TRUE)
              poisson_data_N2_update_clever_full[is.na(poisson_data_N2_update_clever_full)]<-0
              poisson_data_N2_update_clever_full<-poisson_data_N2_update_clever_full%>%select(-H1)
              poisson_data_N2_update_clever_full<-poisson_data_N2_update_clever_full%>%rename(!! paste0("H2.",i):=H2)
            }



  #Evaluate at different point seperately
  canonical_gradient_updated_list<-list()

  for(i in 1:length(h_info_old_list)){

          h_info_old<-h_info_old_list[[i]]

          #Clever covariates info
          obs_id_vec<-lapply(h_info_old, `[[`, 1)
          obs_id_vec<-as.vector(unlist(obs_id_vec))

          time_point_vec<-lapply(h_info_old, `[[`, 2)
          time_point_vec<-as.vector(unlist(time_point_vec))

          H1_vec<-lapply(h_info_old, `[[`, 3)
          H1_vec<-as.vector(unlist(H1_vec))

          H2_vec<-lapply(h_info_old, `[[`, 4)
          H2_vec<-as.vector(unlist(H2_vec))

          clever_covariate_table<-data.frame(obs=obs_id_vec,t=time_point_vec,H1=H1_vec,H2=H2_vec)

          #add clever covariate info to poisson data
          #Poisson data N1
          poisson_data_N1_update_clever <- merge(x=poisson_data_N1_update,y=clever_covariate_table,
                                                 by=c("obs","t"),
                                                 all.x=TRUE)
          poisson_data_N1_update_clever[is.na(poisson_data_N1_update_clever)]<-0
          poisson_data_N1_update_clever<-poisson_data_N1_update_clever%>%select(-H2)

          #Poisson data N2
          poisson_data_N2_update_clever <- merge(x=poisson_data_N2_update,y=clever_covariate_table,
                                                 by=c("obs","t"),
                                                 all.x=TRUE)
          poisson_data_N2_update_clever[is.na(poisson_data_N2_update_clever)]<-0
          poisson_data_N2_update_clever<-poisson_data_N2_update_clever%>%select(-H1)


          #canonical gradient from N1
          poisson_data_N1_update_clever_aug<-cbind(poisson_data_N1_update_clever,part1=poisson_data_N1_update_clever$H1*poisson_data_N1_update_clever$jump)
          poisson_data_N1_update_clever_aug<-cbind(poisson_data_N1_update_clever_aug,part2=
                                                     (1-poisson_data_N1_update_clever$jump)*
                                                     poisson_data_N1_update_clever$timeIntervalLength*
                                                     poisson_data_N1_update_clever$H1*
                                                     exp(as.matrix(cbind(1,poisson_data_N1_update_clever_full[,-c(1:9)]))%*%(as.vector(coef_N1_initial_update))))
          poisson_data_N1_update_clever_aug<-poisson_data_N1_update_clever_aug%>%group_by(obs)%>%summarise(part1=sum(part1),part2=sum(part2))
          poisson_data_N1_update_clever_aug<-poisson_data_N1_update_clever_aug%>%mutate(canonical_N1=part1-part2)

          #canonical gradient from N2
          poisson_data_N2_update_clever_aug<-cbind(poisson_data_N2_update_clever,part1=poisson_data_N2_update_clever$H2*poisson_data_N2_update_clever$jump)
          poisson_data_N2_update_clever_aug<-cbind(poisson_data_N2_update_clever_aug,part2=
                                                     (1-poisson_data_N2_update_clever$jump)*
                                                     poisson_data_N2_update_clever$timeIntervalLength*
                                                     poisson_data_N2_update_clever$H2*
                                                     exp(as.matrix(cbind(1,poisson_data_N2_update_clever_full[,-c(1:9)]))%*%(as.vector(coef_N2_initial_update))))
          poisson_data_N2_update_clever_aug<-poisson_data_N2_update_clever_aug%>%group_by(obs)%>%summarise(part1=sum(part1),part2=sum(part2))
          poisson_data_N2_update_clever_aug<-poisson_data_N2_update_clever_aug%>%mutate(canonical_N2=part1-part2)


          canonical_gradient_updated<-poisson_data_N1_update_clever_aug$canonical_N1+poisson_data_N2_update_clever_aug$canonical_N2
          canonical_gradient_updated_list[[i]]<-canonical_gradient_updated

  }


  return(canonical_gradient_updated_list)
}

