
##Clever Covariates Evaluation
#s1,s2,u are the point of interest for clever covariates
#tt1,tt2,delt1,delt2 are observation i

clever_covaraite_hazards_func_s1_smaller_than_s2<-function(s1,s2,u,tt1,tt2,delt1,delt2,
                                                   A1_time_grid,basis_list_A1_select,coef_A1,
                                                   A2_time_grid,basis_list_A2_select,coef_A2,
                                                   N1_time_grid,basis_list_N1_select,coef_N1,
                                                   N2_time_grid,basis_list_N2_select,coef_N2,
                                                   reverse=F,
                                                   censoring,
                                                   grid_scale_clever_covariate){
      H_1<-0
      H_2<-0

      #Since we also need to do the integration, we need to get the hazard at u for those H that is not 0.
      hazard_N1<-0
      hazard_N2<-0

      ###First case of s1<u<=s2
      if(s1<u & u<=s2){
            ##E21
            E_21<-0

            ##E20
            E_20<-0
            #First scenerio for E_20 to be not 0
            if(tt2>=u & tt1>s1 & tt1<u & delt1==1){
                if(censoring){
                  #prob_denom_E_20<-1
                  # prob_denom_E_20<-prob_denom_E_20*product_integral_func(start=0,end=tt1,0,0,0,A1_time_grid,basis_list_A1_select,coef_A1,type=ifelse(!reverse,"A1","A2"))
                  # prob_denom_E_20<-prob_denom_E_20*product_integral_func(start=0,end=tt1,0,0,0,A2_time_grid,basis_list_A2_select,coef_A2,type=ifelse(!reverse,"A2","A1"))
                  # prob_denom_E_20<-prob_denom_E_20*product_integral_func(start=tt1,end=u,1,1,tt1,A2_time_grid,basis_list_A2_select,coef_A2,type=ifelse(!reverse,"A2","A1"))
                  prob_denom_E_20<-1 #(4-tt1)/4 * (4-u)/4
                }else{
                  prob_denom_E_20<-1
                }


                prob_nume_E_20<-1
                prob_nume_E_20<-prob_nume_E_20*product_integral_func(start=u,end=s2,1,1,tt1,N2_time_grid,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))


                E_20<-prob_nume_E_20/prob_denom_E_20


                hazard_N1<-0
                hazard_N2<-hazard_func(u,1,1,tt1,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))


            }

            #Second scenerio for E_20 to be not 0 (NOTE THAT WE MOVE THE EDGE CASE u=tt1 to the second scenerio to make the E10 calculation trival)
            if(tt2>=u & tt1>=u){
                if(censoring){
                  # prob_denom_E_20<-1
                  # prob_denom_E_20<-prob_denom_E_20*product_integral_func(start=0,end=u,0,0,0,A1_time_grid,basis_list_A1_select,coef_A1,type=ifelse(!reverse,"A1","A2"))
                  # prob_denom_E_20<-prob_denom_E_20*product_integral_func(start=0,end=u,0,0,0,A2_time_grid,basis_list_A2_select,coef_A2,type=ifelse(!reverse,"A2","A1"))
                  prob_denom_E_20<-1 #(4-u)/4 * (4-u)/4
                }else{
                  prob_denom_E_20<-1
                }

                prob_nume_E_20_part2<-1
                prob_nume_E_20_part2<-prob_nume_E_20_part2*product_integral_func(start=u,end=s2,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type=ifelse(!reverse,"N1","N2"))
                prob_nume_E_20_part2<-prob_nume_E_20_part2*product_integral_func(start=u,end=s2,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))


                prob_nume_E_20_part1<-0
                if(u<s2){
                  seq_temp<-sort(unique(c(seq(u,s2,by=grid_scale_clever_covariate),s2)))
                  left_point_eval<-sapply(seq_temp[-length(seq_temp)],function(x){
                        hazard_func(x,0,0,0,basis_list_N1_select,coef_N1,type=ifelse(!reverse,"N1","N2"))*
                        product_integral_func(start=u,end=x,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type=ifelse(!reverse,"N1","N2"))*
                        product_integral_func(start=u,end=x,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))*
                        product_integral_func(start=x,end=s2,1,1,x,N2_time_grid,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))
                        })
                  length_eval<-diff(seq_temp)
                  prob_nume_E_20_part1<-sum(left_point_eval*length_eval)
                }

                E_20<-(prob_nume_E_20_part1+prob_nume_E_20_part2)/prob_denom_E_20

                hazard_N1<-hazard_func(u,0,0,0,basis_list_N1_select,coef_N1,type=ifelse(!reverse,"N1","N2"))
                hazard_N2<-hazard_func(u,0,0,0,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))
            }

            ##E10
            E_10<-0
            if(tt2>=u & tt1>=u){
              E_10<-E_20
            }

            #E11
            E_11<-0
            if(tt2>=u & tt1>=u){
              prob_denom_E_11<-prob_denom_E_20

              prob_nume_E_11<-1
              prob_nume_E_11<-prob_nume_E_11*product_integral_func(start=u,end=s2,1,1,u,N2_time_grid,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))

              E_11<-prob_nume_E_11/prob_denom_E_11
            }

            #clever covariates
            H_1<-E_11-E_10
            H_2<-E_21-E_20
      }

      ###Second case of u<=s1<=s2
      if(u<=s1 & s1<=s2){
        ##E21
        E_21<-0

        ##E20
        E_20<-0
        if(tt2>=u & tt1>=u){
          if(censoring){
            # prob_denom_E_20<-1
            # prob_denom_E_20<-prob_denom_E_20*product_integral_func(start=0,end=u,0,0,0,A1_time_grid,basis_list_A1_select,coef_A1,type=ifelse(!reverse,"A1","A2"))
            # prob_denom_E_20<-prob_denom_E_20*product_integral_func(start=0,end=u,0,0,0,A2_time_grid,basis_list_A2_select,coef_A2,type=ifelse(!reverse,"A2","A1"))
            prob_denom_E_20<-1 #(4-u)/4 * (4-u)/4
          }else{
            prob_denom_E_20<-1
          }

          prob_nume_E_20_part2<-1
          prob_nume_E_20_part2<-prob_nume_E_20_part2*product_integral_func(start=u,end=s2,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type=ifelse(!reverse,"N1","N2"))
          prob_nume_E_20_part2<-prob_nume_E_20_part2*product_integral_func(start=u,end=s2,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))


          prob_nume_E_20_part1<-0
          if(s1<s2){
            seq_temp<-sort(unique(c(seq(s1,s2,by=grid_scale_clever_covariate),s2)))
            left_point_eval<-sapply(seq_temp[-length(seq_temp)],function(x){
              hazard_func(x,0,0,0,basis_list_N1_select,coef_N1,type=ifelse(!reverse,"N1","N2"))*
                product_integral_func(start=u,end=x,0,0,0,N1_time_grid,basis_list_N1_select,coef_N1,type=ifelse(!reverse,"N1","N2"))*
                product_integral_func(start=u,end=x,0,0,0,N2_time_grid,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))*
                product_integral_func(start=x,end=s2,1,1,x,N2_time_grid,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))
            })
            length_eval<-diff(seq_temp)
            prob_nume_E_20_part1<-sum(left_point_eval*length_eval)
          }

          E_20<-(prob_nume_E_20_part1+prob_nume_E_20_part2)/prob_denom_E_20

          hazard_N1<-hazard_func(u,0,0,0,basis_list_N1_select,coef_N1,type=ifelse(!reverse,"N1","N2"))
          hazard_N2<-hazard_func(u,0,0,0,basis_list_N2_select,coef_N2,type=ifelse(!reverse,"N2","N1"))
        }

        ##E10
        E_10<-0
        if(tt2>=u & tt1>=u){
          E_10<-E_20
        }

        #E11
        E_11<-0


        #clever covariates
        H_1<-E_11-E_10
        H_2<-E_21-E_20
      }

      if(!reverse){
        H_1<-H_1
        H_2<-H_2
        hazard_N1<-hazard_N1
        hazard_N2<-hazard_N2
      }else{
        H_1<-H_2
        H_2<-H_1
        hazard_N1<-hazard_N2
        hazard_N2<-hazard_N1
      }
      return(c(H_1,H_2,hazard_N1,hazard_N2))
}




clever_covaraite_hazards_func<-function(s1,s2,u,tt1,tt2,delt1,delt2,
                                  A1_time_grid,basis_list_A1_select,coef_A1,
                                  A2_time_grid,basis_list_A2_select,coef_A2,
                                  N1_time_grid,basis_list_N1_select,coef_N1,
                                  N2_time_grid,basis_list_N2_select,coef_N2,
                                  censoring,
                                  grid_scale_clever_covariate
                                 ){

    if(s1<=s2){
      result<-clever_covaraite_hazards_func_s1_smaller_than_s2(s1,s2,u,tt1,tt2,delt1,delt2,
                                                                A1_time_grid,basis_list_A1_select,coef_A1,
                                                                A2_time_grid,basis_list_A2_select,coef_A2,
                                                                N1_time_grid,basis_list_N1_select,coef_N1,
                                                                N2_time_grid,basis_list_N2_select,coef_N2,
                                                                reverse=F,
                                                               censoring=censoring,
                                                               grid_scale_clever_covariate)
    }else{
      result<-clever_covaraite_hazards_func_s1_smaller_than_s2(s2,s1,u,tt2,tt1,delt2,delt1,
                                                               A2_time_grid,basis_list_A2_select,coef_A2,
                                                               A1_time_grid,basis_list_A1_select,coef_A1,
                                                               N2_time_grid,basis_list_N2_select,coef_N2,
                                                               N1_time_grid,basis_list_N1_select,coef_N1,
                                                               reverse=T,
                                                               censoring=censoring,
                                                               grid_scale_clever_covariate)
    }
   return(result)
}



