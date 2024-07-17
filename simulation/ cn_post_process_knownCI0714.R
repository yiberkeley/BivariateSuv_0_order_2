

library(ggplot2)
library(ggpubr)
library(purrr)
library(dplyr)
library(plotly)
library(foreach)
library(doParallel)
library(future.apply)
library(fMultivar)
source("sim_data.R")
source("sim_data_2.R")
`%+%` <- function(a, b) paste0(a, b)

#START

#Input the simulation result
input_result_fun<-function(filename_first,seed_vals){
  result<-list()
  for(i in 1:length(seed_vals)){
    result_temp<-readRDS(filename_first %+% seed_vals[i] %+% ".RDS")
    result[[i]]<-result_temp
  }
  return(result)
}



result.l3.targeted.1 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
                                        seed_vals =  c(101,500,2000,4321,10113,45000,48001,50002))


#knots, SCS, up end, trunc
result.l3.targeted.1.knotTrucUpEnd.1 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002))
result.l3.targeted.1.knotTrucUpEnd.2 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+1111)
result.l3.targeted.1.knotTrucUpEnd.3 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+2222)
result.l3.targeted.1.knotTrucUpEnd.4 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+3333)
# result.l3.targeted.1.knotTrucUpEnd.5 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+4444)
# result.l3.targeted.1.knotTrucUpEnd.6 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+5555)
# result.l3.targeted.1.knotTrucUpEnd.7 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+6666)
# result.l3.targeted.1.knotTrucUpEnd.8 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+7777)
# result.l3.targeted.1.knotTrucUpEnd.9 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_200_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+8888)

result.l3.targeted.1.knotTrucUpEnd<-c(result.l3.targeted.1.knotTrucUpEnd.1,
                                      result.l3.targeted.1.knotTrucUpEnd.2,
                                      result.l3.targeted.1.knotTrucUpEnd.3,
                                      result.l3.targeted.1.knotTrucUpEnd.4)
                                      # result.l3.targeted.1.knotTrucUpEnd.5,
                                      # result.l3.targeted.1.knotTrucUpEnd.6,
                                      # result.l3.targeted.1.knotTrucUpEnd.7,
                                      # result.l3.targeted.1.knotTrucUpEnd.8,
                                      # result.l3.targeted.1.knotTrucUpEnd.9)

#knots, SCS, up end, trunc
result.l3.targeted.2.knotTrucUpEnd.1 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_800_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002))
result.l3.targeted.2.knotTrucUpEnd.2 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_800_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+1111)
result.l3.targeted.2.knotTrucUpEnd.3 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_800_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+2222)
result.l3.targeted.2.knotTrucUpEnd.4 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_800_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+3333)
# result.l3.targeted.2.knotTrucUpEnd.5 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_800_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+4444)
# result.l3.targeted.2.knotTrucUpEnd.6 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_800_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+5555)
# result.l3.targeted.2.knotTrucUpEnd.7 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_800_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+6666)
# result.l3.targeted.2.knotTrucUpEnd.8 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_800_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+7777)
# result.l3.targeted.2.knotTrucUpEnd.9 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_800_",
#                                                         seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+8888)

result.l3.targeted.2.knotTrucUpEnd<-c(result.l3.targeted.2.knotTrucUpEnd.1,
                                      result.l3.targeted.2.knotTrucUpEnd.2,
                                      result.l3.targeted.2.knotTrucUpEnd.3,
                                      result.l3.targeted.2.knotTrucUpEnd.4)
                                      # result.l3.targeted.2.knotTrucUpEnd.5,
                                      # result.l3.targeted.2.knotTrucUpEnd.6,
                                      # result.l3.targeted.2.knotTrucUpEnd.7,
                                      # result.l3.targeted.2.knotTrucUpEnd.8,
                                      # result.l3.targeted.2.knotTrucUpEnd.9)



result.l3.targeted.4.knotTrucUpEnd.1 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_improve1800_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002))
result.l3.targeted.4.knotTrucUpEnd.2 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_improve1800_",
                                                        seed_vals = c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+1111)
result.l3.targeted.4.knotTrucUpEnd.3 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_improve1800_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+2222)
result.l3.targeted.4.knotTrucUpEnd.4 <-input_result_fun(filename_first="071424_various/knT_071424_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_improve1800_",
                                                        seed_vals =  c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)+3333)

result.l3.targeted.4.knotTrucUpEnd<-c(result.l3.targeted.4.knotTrucUpEnd.1,
                                      result.l3.targeted.4.knotTrucUpEnd.2,
                                      result.l3.targeted.4.knotTrucUpEnd.3,
                                      result.l3.targeted.4.knotTrucUpEnd.4)
t1<-result.l3.targeted.1[[1]][[1]][[1]]$t1
t2<-result.l3.targeted.1[[1]][[1]][[1]]$t2


#Evaluate the truth
# t1_plot<-t1
# t2_plot<-t2
# index<-seq(1,length(t1)*length(t2))

t1_plot<-t1[seq(1,10000,by=500)]
t2_plot<-t2[seq(1,10000,by=500)]
index<-which(t1%in%t1_plot & t2%in%t2_plot)

datgen<-function(n){return(datgen_knownIntensity_no_censoring(n,1000))}
plot_data<-expand.grid(t1=t1_plot,t2=t2_plot)
set.seed(100)
nCores<-detectCores()-6
plan(multisession,workers = nCores)
data_full_sim<-datgen(20000)

truth.fineGridGen<-future_apply(plot_data,1,function(x) mean(data_full_sim$t1<=x[1] & data_full_sim$t2<=x[2]))

temp_dataTrue<-plot_data%>%mutate(cdf=truth.fineGridGen)
temp_dataTrue<-temp_dataTrue%>%arrange(t1,t2)
####
#Simulate the data and get the empirical histogram

#Investigation into the fit

plot_data_func<-function(result1,t1,t2,truth1,index,update){
  result_temp<-result1
  sum_cdf_initial<-  result_temp[[1]][[1]][[7]]$cdf_curve[index]
  if(update){sum_cdf_update<-   result_temp[[1]][[2]][[7]]$cdf_curve[index]}
  for(i in 2: length(result_temp)){
    sum_cdf_initial<-sum_cdf_initial+result_temp[[i]][[1]][[7]]$cdf_curve[index]
    if(update){sum_cdf_update<-sum_cdf_update+result_temp[[i]][[2]][[7]]$cdf_curve[index]}
  }
  mean_cdf_initial1<-sum_cdf_initial/length(result_temp)
  if(update){mean_cdf_update1<-sum_cdf_update/length(result_temp)}

  rss_cdf_initial<-(result_temp[[1]][[1]][[7]]$cdf_curve[index]-mean_cdf_initial1)^2
  if(update){rss_cdf_update<-(result_temp[[1]][[2]][[7]]$cdf_curve[index]-mean_cdf_update1)^2}
  for(i in 2: length(result_temp)){
    rss_cdf_initial<-rss_cdf_initial+(result_temp[[i]][[1]][[7]]$cdf_curve[index]-mean_cdf_initial1)^2
    if(update){rss_cdf_update<-rss_cdf_update+(result_temp[[i]][[2]][[7]]$cdf_curve[index]-mean_cdf_update1)^2}
  }
  sd_cdf_initial1<-sqrt(rss_cdf_initial/length(result_temp))
  if(update){sd_cdf_update1<-sqrt(rss_cdf_update/length(result_temp))}



  plot_data1.1<-data_frame(t1=t1[index],t2=t2[index],mean_cdf_bias=mean_cdf_initial1, sd_surve=sd_cdf_initial1,update=0)
  if(update){plot_data1.2<-data_frame(t1=t1[index],t2=t2[index],mean_cdf_bias=mean_cdf_update1 , sd_surve=sd_cdf_update1,update=1)}
  if(update){
    return(rbind(plot_data1.1,plot_data1.2))
  }else{
    return(rbind(plot_data1.1))
  }
}


#########

result1<-result.l3.targeted.4.knotTrucUpEnd
#result1<-result.l3.targeted.1.knotFineGridTruc
truth1<-truth.fineGridGen
update<-F

result_temp1<-result1
result_temp<-NULL
for(i in 1:length(result_temp1)){
  if(result_temp1[[i]][[8]]==result_temp1[[i]][[9]] &result_temp1[[i]][[8]]=="optimal"){
    result_temp<-c(result_temp,i)
  }
}
result1<-result1[result_temp]

library(plotly)
plot_data<-plot_data_func(result1,t1,t2,truth1,index,update)
plot_data<-plot_data%>%arrange(t1,t2)
plot_data$mean_cdf_bias<-plot_data$mean_cdf_bias-temp_dataTrue$cdf

plot_ly(plot_data,x=~t1,y=~t2,z=~mean_cdf_bias,color=~update,colors = c('#BF382A', '#0C4B8E'))
plot_ly(plot_data,x=~t1,y=~t2,z=~temp_dataTrue$cdf,color=~update,colors = c('#BF382A', '#0C4B8E'))
plot_ly(plot_data,x=~t1,y=~t2,z=~mean_cdf_bias+temp_dataTrue$cdf,color=~update,colors = c('#BF382A', '#0C4B8E'))
plot_ly(plot_data,x=~t1,y=~t2,z=~mean_cdf_bias,color=~update,colors = c('#BF382A', '#0C4B8E')) #Good region

plot_ly(plot_data%>%filter(t1<=0.701 & t2<=0.701),x=~t1,y=~t2,z=~mean_cdf_bias,color=~update,colors = c('#BF382A', '#0C4B8E')) #Good region
plot_ly(plot_data%>%filter(t1<=0.701 & t2<=0.701),x=~t1,y=~t2,z=~abs(mean_cdf_bias),color=~update,colors = c('#BF382A', '#0C4B8E')) #Good region
plot_ly(plot_data%>%filter(t1<=0.701 & t2<=0.701),x=~t1,y=~t2,z=~sd_surve,color=~update,colors = c('#BF382A', '#0C4B8E')) #Good region

#bias/sd ratio
plot_ly(plot_data%>%filter(t1<=0.7 & t2<=0.7, update==0),x=~t1,y=~t2,z=~abs(mean_cdf_bias)/sd_surve,color=~update,colors = c('#BF382A', '#0C4B8E')) #Good region
plot_ly(plot_data%>%filter(update==0),x=~t1,y=~t2,z=~abs(mean_cdf_bias)/sd_surve,color=~update,colors = c('#BF382A', '#0C4B8E'))

plot_dataGR<-plot_data%>%filter(t1<=0.7 & t2<=0.7, update==0)
hist(abs(plot_dataGR$mean_cdf_bias)/plot_dataGR$sd_surve,breaks=100,main="MLE, parametric DGP")

mean(abs(plot_data$mean_cdf_bias)/plot_data$sd_surve<=1,na.rm=T)

###########################################
###########################################
#N1

#result_temp1<-result.l3.targeted.1.knotTrucUnifFineGrid
result_temp1<-result.l3.targeted.4.knotTrucUpEnd
result_temp<-NULL
for(i in 1:length(result_temp1)){
  if(1==1){#result_temp1[[i]][[8]]==result_temp1[[i]][[9]] &result_temp1[[i]][[8]]=="optimal"){
    temp_data<-data_frame(coef=c(1:length(coef_t1)),value=result_temp1[[i]][[2]],round=i)
    result_temp<-rbind(result_temp,temp_data)
  }
}

sum_optimal<-0
for(i in 1:length(result_temp1)){
  if(result_temp1[[i]][[8]]=="optimal"){sum_optimal=sum_optimal+1}
}
sum_optimal
length(result_temp1)

result_true<-data_frame(coef=c(1:length(coef_t1)),value=coef_t1,round='true')
ggplot(result_temp,aes(x=coef,y=value,group=round))+
  geom_point(aes(color=round))+
  geom_line(aes(color=round))+
  geom_point(data = result_true, col = 'red')+
  ggtitle("N1_coef")

ggplot(result_temp,aes(x=coef,y=value,group=round))+
  geom_point(aes(color=round))+
  geom_line(aes(color=round))+
  geom_point(data = result_true, col = 'red')+
  ggtitle("N1_coef")+ylim(c(-4,4))

ggplot(result_temp,aes(x=value,group=round))+
  geom_histogram(aes(color=round))+
  facet_wrap(.~coef)+
  geom_vline(data = result_true, aes(xintercept=value),col = 'red')+
  ggtitle("N1_coef")

summary_table<-result_temp%>%group_by(coef)%>%summarise(mean=mean(value),sd=sd(value))
summary_table$mean-coef_t1
summary_table$sd
abs(summary_table$mean-coef_t1)/summary_table$sd

beta_estimates<-(result_temp%>%filter(coef==1))$value
qqnorm(beta_estimates, pch = 1, frame = FALSE)
qqline(beta_estimates, col = "steelblue", lwd = 2)

beta_estimates<-(result_temp%>%filter(coef==2))$value
qqnorm(beta_estimates, pch = 1, frame = FALSE)
qqline(beta_estimates, col = "steelblue", lwd = 2)


#N2
result_temp1<-result.l3.targeted.4.knotTrucUpEnd
result_temp<-NULL
for(i in 1:length(result_temp1)){
  if(1==1){#result_temp1[[i]][[8]]==result_temp1[[i]][[9]] &result_temp1[[i]][[8]]=="optimal"){
    temp_data<-data_frame(coef=c(1:length(coef_t2)),value=result_temp1[[i]][[4]],round=i)
    result_temp<-rbind(result_temp,temp_data)
  }
}

sum_optimal<-0
for(i in 1:length(result_temp1)){
  if(result_temp1[[i]][[9]]=="optimal"){sum_optimal=sum_optimal+1}
}
sum_optimal
length(result_temp1)

result_true<-data_frame(coef=c(1:length(coef_t2)),value=coef_t2,round='true')
ggplot(result_temp,aes(x=coef,y=value,group=round))+
  geom_point(aes(color=round))+
  geom_line(aes(color=round))+
  geom_point(data = result_true, col = 'orange')+
  ggtitle("N2_coef")

ggplot(result_temp,aes(x=value,group=round))+
  geom_histogram(aes(color=round))+
  facet_wrap(.~coef)+
  geom_vline(data = result_true, aes(xintercept=value),col = 'red')+
  ggtitle("N2_coef")

summary_table<-result_temp%>%group_by(coef)%>%summarise(mean=median(value),sd=sd(value))
summary_table$mean-coef_t2
summary_table$sd
abs(summary_table$mean-coef_t2)/summary_table$sd

beta_estimates<-(result_temp%>%filter(coef==1))$value
qqnorm(beta_estimates, pch = 1, frame = FALSE)
qqline(beta_estimates, col = "steelblue", lwd = 2)

beta_estimates<-(result_temp%>%filter(coef==2))$value
qqnorm(beta_estimates, pch = 1, frame = FALSE)
qqline(beta_estimates, col = "steelblue", lwd = 2)

##Looks like ECOS better than the SCS
##Looks like FineGrid and NoFineGrid roughly same
##Looks like with knots is slightly better than no knots

#Trucation seems to introduce the downward trend, but parameters estimation is better



##To do

##Step1:Remove bias in the betas
#CHOICE 1 (First try)
#choose the points uniformly in the small intervals.(Done, no huge improvements)

#CHOICE 2
#repeated data at the same discretization with end points.

#Discretized truth.


##Step2:Make sure the density is exactly evaluated at the fitted betas (0 order case)


##
##Step3:
#1.Maybe very fine grid.
#2.Or (PRIORITY)"piecewise linear approximation through a 1st order hal(unpenalized) fit" (estimated density fitted on the tensor products of the 1 order indicators to get the approximations.)
#A cts function goes through all the estimated density points of a grid.
#Bias level is 1/h^2, so even when h=20, we are better off than 0.0025
#(t1-u1)*1(t1>u1) (t2-u2)*1(t2>u2)

#3.(PRIORITY)Calculated the fitted cdf explicitly from the fitted betas, no grid approximation needed.(Oracle benchmark)

