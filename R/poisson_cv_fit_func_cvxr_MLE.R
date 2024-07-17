#This is a MLE version for investigaation

#Poisson cv fit helper function
poisson_cv_fit_func_cvxr<-function(poisson_data_N3,basis_list_N3,nfolds,penalty.factor_N3,undersmooth=F,l1_norm_start=F,uniform_ratio=F,
                                   offset_initial=NULL,length_clever=NULL,coef_initial=NULL,offset_relax=NULL){

    offset<-log(poisson_data_N3$timeIntervalLength)
    weights<- poisson_data_N3$weights
    y<-poisson_data_N3$jump
    x_basis_poisson<-as.matrix(poisson_data_N3[,-c(1:9)])
    time_vec<-as.vector(poisson_data_N3$t)

    if(is.null(offset_initial)){
      fit_N3<-cvrx_poisson_l1(1000000,x_basis_poisson,offset,weights,y,penalty.factor_N3,time_vec,uniform_ratio)
    }else{
      fit_N3<-cvrx_poisson_l1_offset(1000000,x_basis_poisson,offset,weights,y,penalty.factor_N3,time_vec,uniform_ratio,length_clever,coef_initial)
    }


  beta<-fit_N3[[1]][[1]]
  beta_0<-fit_N3[[1]][[2]]
  cv_selected_fit_risk_N3<-fit_N3[[1]][[3]]
  solver_status<-fit_N3[[1]][[4]]

  coef_N3_initial_list<-list()
  basis_list_N3_select_list<-list()
  penalty_N3_select_list<-list()

  #beta<-ifelse(beta<=10^-4,0,beta) #reduce the beta to 0 if solved to be very close to 0

  #Return as 0 betas as 0
  for(i in 1:1){
    coef_N3_initial_list[[i]]<-c(beta_0,beta)#[abs(beta)>10^-4])
    basis_list_N3_select_list[[i]]<-basis_list_N3#[abs(beta)>10^-4]
    penalty_N3_select_list[[i]]<-penalty.factor_N3#[abs(beta)>10^-4]
  }



  result_list<-list()
  result_list[[1]]<-1
  result_list[[2]]<-1
  result_list[[3]]<-cv_selected_fit_risk_N3
  result_list[[4]]<-coef_N3_initial_list
  result_list[[5]]<-basis_list_N3_select_list
  result_list[[6]]<-penalty_N3_select_list
  result_list[[7]]<-solver_status

  return(result_list)
}
