print("data_3")
#Specify the known conditional intensity when there is no censoring
#The conditional intensity specification for N1
set.seed(1000)
library(Rlab)
t1_gen<-runif(20,min=0,max=1)
t2_happened_gen<-rbern(20,0.5)
t2_happened_t2_gen<-t2_happened_gen*runif(20,min=0,max=1)*t1_gen
t2_happened_delta_gen<-t2_happened_gen*rbern(20,0.7)
data_t1_basis_gen<-data.frame(t=t1_gen,t2_happened=as.integer(t2_happened_gen),t2_happened_delta=as.integer(t2_happened_delta_gen),t2_happened_t2=t2_happened_t2_gen)
#Specify the basis degree and the num_knots.
#In this case the hazard are piecewise linear in each strata and cts.
basis_list_t1 <- enumerate_basis(data_t1_basis_gen,
                                 max_degree = 2,
                                 smoothness_orders = c(0,0,0,0),#c(1,0,0,1),
                                 #include_zero_order  = T,
                                 #include_lower_order = T,
                                 num_knots = c(3,NULL,NULL,3)
)
#Process the basis to get rid of the 0 points which is really reduced to lower interactions and the intercept in the 0 interactions.
basis_list_t1_process<-NULL
j=1
for(i in 1:length(basis_list_t1)){
  if(sum(basis_list_t1[[i]]$cutoffs==0)==0){
    basis_list_t1_process[[j]]<-basis_list_t1[[i]]
    j=j+1
  }
}
basis_list_t1<-basis_list_t1_process
#Process again to get rid of non-sense interaction t2_happened and t2_happened_delta
basis_list_t1_process<-NULL
j=1
for(i in 1:length(basis_list_t1)){
  if(sum(basis_list_t1[[i]]$cols==c(2,3))!=2){
    basis_list_t1_process[[j]]<-basis_list_t1[[i]]
    j=j+1
  }
}
basis_list_t1<-basis_list_t1_process
#Process again to get rid of non-sense interaction t2_happened and t2_happened_tt2
basis_list_t1_process<-NULL
j=1
for(i in 1:length(basis_list_t1)){
  if(sum(basis_list_t1[[i]]$cols==c(2,4))!=2){
    basis_list_t1_process[[j]]<-basis_list_t1[[i]]
    j=j+1
  }
}
basis_list_t1<-basis_list_t1_process
basis_list_t1<-basis_list_t1[c(2:5,7,9,11,16)]
#Process to change the extreme values near 0 or 1 to increase the identifiability
basis_list_t1[[7]]$cutoffs<-c(0.2,1)
#Process to get rid of the type info since it is always uncensored in our simulation.
basis_list_t1<-basis_list_t1[-c(4,7)]
#coef_t1<-c(1.4,1.76,0.2,0.15,0.5,1.2,0.45)
coef_t1<-runif(length(basis_list_t1)+1,0.3,1)
coef_t1[1]<-0.88
coef_t1[3]<-1.3
coef_t1[4]<-0.38
coef_t1[5]<0.9
coef_t1[6]<-0.5

################
##For investigation purpose, we will test out several simple cases
#No interaction btw jumping processes and only one knot points.
basis_list_t1[[1]]$cutoffs<-c(0.2)
basis_list_t1[[2]]$cutoffs<-c(0.33)
basis_list_t1[[4]]$cutoffs<-c(0.11)
basis_list_t1[[5]]$cutoffs<-c(0.25,1)
basis_list_t1[[6]]$cutoffs<-c(0.15,0.35)


################

#The conditional intensity specification for N2
t2_gen<-runif(20,min=0,max=1)
t1_happened_gen<-rbern(20,0.3)
t1_happened_t1_gen<-t1_happened_gen*runif(20,min=0,max=1)*t2_gen
t1_happened_delta_gen<-t1_happened_gen*rbern(20,0.7)
data_t2_basis_gen<-data.frame(t=t2_gen,t1_happened=t1_happened_gen,t1_happened_delta=t1_happened_delta_gen,t1_happened_t1=t1_happened_t1_gen)
basis_list_t2 <- enumerate_basis(data_t2_basis_gen,
                                 max_degree = 2,
                                 smoothness_orders = c(0,0,0,0),#c(1,0,0,1),
                                 #include_zero_order  = T,
                                 #include_lower_order = T,
                                 num_knots = c(3,NULL,NULL,3)
)
#Process the basis to get rid of the 0 points which is really reduced to lower interactions and the intercept in the 0 interactions.
basis_list_t2_process<-NULL
j=1
for(i in 1:length(basis_list_t2)){
  if(sum(basis_list_t2[[i]]$cutoffs==0)==0){
    basis_list_t2_process[[j]]<-basis_list_t2[[i]]
    j=j+1
  }
}
basis_list_t2<-basis_list_t2_process
#Process again to get rid of non-sense interaction t2_happened and t2_happened_delta
basis_list_t2_process<-NULL
j=1
for(i in 1:length(basis_list_t2)){
  if(sum(basis_list_t2[[i]]$cols==c(2,3))!=2){
    basis_list_t2_process[[j]]<-basis_list_t2[[i]]
    j=j+1
  }
}
basis_list_t2<-basis_list_t2_process
#Process again to get rid of non-sense interaction t2_happened and t2_happened_tt2
basis_list_t2_process<-NULL
j=1
for(i in 1:length(basis_list_t2)){
  if(sum(basis_list_t2[[i]]$cols==c(2,4))!=2){
    basis_list_t2_process[[j]]<-basis_list_t2[[i]]
    j=j+1
  }
}
basis_list_t2<-basis_list_t2_process
basis_list_t2<-basis_list_t2[c(2,4,5,6,8)]
#Process to get rid of the type info since it is always uncensored in our simulation.
basis_list_t2<-basis_list_t2[-c(3)]
coef_t2<-runif(length(basis_list_t2)+1,0.3,1)
coef_t2[1]<-1.2
coef_t2[2]<-0.3
coef_t2[4]<-0.8

################
##For investigation purpose, we will test out several simple cases
#No interaction btw jumping processes and only one knot points.
basis_list_t2[[1]]$cutoffs<-c(0.17)
basis_list_t2[[3]]$cutoffs<-c(0.3)
basis_list_t2[[4]]$cutoffs<-c(0.35,1)
################



gen_t1_t2<-function(c1,c2){
  grid_gen<-seq(0,1,by=0.001) #0.001 to the fineGridGen.
  #Make sure the landmarks are just 100 maybe

  t1_happened=0
  t1_happened_t1=0
  t1_happened_delta=0
  t1_full=0
  t1_full_t1=0

  t2_happened=0
  t2_happened_t2=0
  t2_happened_delta=0
  t2_full=0
  t2_full_t2=0
  for(i in 1:(length(grid_gen)-1)){

    if(t1_full==0){
      #Jumping process N1
      #Use midpoint for the prob
      #temp<-data.frame(t=(grid_gen[i]+grid_gen[i+1])/2,t2_happened=t2_happened,t2_happened_delta=t2_happened_delta,t2_happened_t2=t2_happened_t2)
      #Use lower point for the prob
      temp<-data.frame(t=grid_gen[i],t2_happened=t2_happened,t2_happened_delta=t2_happened_delta,t2_happened_t2=t2_happened_t2)

      temp2<-make_design_matrix(as.matrix(temp),basis_list_t1)
      temp2<-c(1,as.vector(temp2))
      log_hazard_midpoint_N1<-sum(coef_t1*temp2)
      prob_interval_N1<-exp(-(grid_gen[i+1]-grid_gen[i])*exp(log_hazard_midpoint_N1)) #Not jumping
      jump_interval_N1<-rbern(1,1-prob_interval_N1)
      if(t1_happened==0){
        if(jump_interval_N1==1){
          t1_happened=1
          t1_happened_t1=grid_gen[i+1]#runif(1,grid_gen[i],grid_gen[i+1])#(grid_gen[i]+grid_gen[i+1])/2
          t1_happened_delta=1
          t1_full=1
          t1_full_t1=t1_happened_t1
        }else if(c1<grid_gen[i+1]){
          t1_happened=1
          t1_happened_t1=c1
          t1_happened_delta=0
        }else{
          t1_happened=0
          t1_happened_t1=0
          t1_happened_delta=0
        }
      }else{
        if(jump_interval_N1==1){
          t1_full=1
          t1_full_t1=grid_gen[i+1]#runif(1,grid_gen[i],grid_gen[i+1])#(grid_gen[i]+grid_gen[i+1])/2
        }
      }
    }

    if(t2_full==0){
      #Jumping process N2
      #temp<-data.frame(t=(grid_gen[i]+grid_gen[i+1])/2,t1_happened=t1_happened,t1_happened_delta=t1_happened_delta,t1_happened_t1=t1_happened_t1)
      #Use lower point for the prob
      temp<-data.frame(t=grid_gen[i],t1_happened=t1_happened,t1_happened_delta=t1_happened_delta,t1_happened_t1=t1_happened_t1)

      temp2<-make_design_matrix(as.matrix(temp),basis_list_t2)
      temp2<-c(1,as.vector(temp2))
      log_hazard_midpoint_N2<-sum(coef_t2*temp2)
      prob_interval_N2<-exp(-(grid_gen[i+1]-grid_gen[i])*exp(log_hazard_midpoint_N2)) #Not jumping
      jump_interval_N2<-rbern(1,1-prob_interval_N2)
      if(t2_happened==0){
        if(jump_interval_N2==1){
          t2_happened=1
          t2_happened_t2=grid_gen[i+1]#runif(1,grid_gen[i],grid_gen[i+1])#(grid_gen[i]+grid_gen[i+1])/2
          t2_happened_delta=1
          t2_full=1
          t2_full_t2=t2_happened_t2
        }else if(c2<grid_gen[i+1]){
          t2_happened=1
          t2_happened_t2=c2
          t2_happened_delta=0
        }else{
          t2_happened=0
          t2_happened_t2=0
          t2_happened_delta=0
        }
      }else{
        if(jump_interval_N2==1){
          t2_full=1
          t2_full_t2=grid_gen[i+1]#runif(1,grid_gen[i],grid_gen[i+1])#(grid_gen[i]+grid_gen[i+1])/2
        }
      }
    }
  }

  if(t1_full==0){
    t1_full_t1=1
  }

  if(t2_full==0){
    t2_full_t2=1
  }

  return(c(t1_full_t1,t2_full_t2))
}




datgen_knownIntensity_no_censoring<-function(n,seed){
  c1<-rep(4,n)
  c2<-rep(4,n)

  c1<-ceiling(c1*1000)/1000
  c2<-ceiling(c1*1000)/1000

  set.seed(seed)
  t1<-rep(0,n)
  t2<-rep(0,n)
  for(i in 1:n){
    t1_t2<-gen_t1_t2(c1[i],c2[i])
    t1[i]<-t1_t2[1]
    t2[i]<-t1_t2[2]
  }

  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)
  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}

datgen_knownIntensity_no_censoring_trunc<-function(n,seed){
  c1<-rep(4,n)
  c2<-rep(4,n)

  c1<-ceiling(c1*1000)/1000
  c2<-ceiling(c1*1000)/1000

  set.seed(seed)
  t1<-rep(0,n)
  t2<-rep(0,n)
  for(i in 1:n){
    t1_t2<-gen_t1_t2(c1[i],c2[i])
    t1[i]<-t1_t2[1]
    t2[i]<-t1_t2[2]
  }

  tt2<-apply(cbind(t2,c2),1,min)
  tt1<- apply(cbind(t1,c1),1,min)
  delt1<-  as.numeric(t1<=c1)
  delt2<- as.numeric(t2<=c2)

  index_temp<-which(tt1>0.68)
  tt1[index_temp]<-0.68
  delt1[index_temp]<-0

  index_temp<-which(tt2>0.68)
  tt2[index_temp]<-0.68
  delt2[index_temp]<-0


  type<-1
  w<-0
  df <- data.frame(t1 = t1, t2 = t2,tt1 = tt1, tt2 = tt2,c2=c2, delt1=delt1,delt2=delt2,w,type)
  return(df)
}

