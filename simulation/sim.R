source("utils.R")
source("utils_MLE.R")

#n=500, perc_sample=0.1,perc_sample_basis=0.1
#n=1000, monitoring times: perc_sample=0.02 & full basis subsampling: perc_sample_basis=0.02
#n=2000, perc_sample=0.01,perc_sample_basis=0.01
#seed_vals<-c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)
seed_vals<-c(101,500,2000,4321,8688)
#seed_vals<-c(10113,36899,45000,48001,50002)+3333


for (seed in seed_vals) {
  n=800
  res <- run_simMLE(n = n,
                 seed = seed,
                 type = "knT",
                 nfolds=3,
                 cv_choice="lambda",
                 undersmooth_type="cv_all",
                 basis_option="0_order",
                 masking_option="no_masking", #0.7
                 weight_option="no_tail_shrinkage",
                 penalty_type="usual",
                 undersmooth=F,
                 targeting=F,
                 censoring=F,
                 relaxed=F,
                 relaxed_initial=F,
                 offset_initial=T,
                 offset_relax=F,
                 perc_sample=1,
                 perc_sample_basis=1,
                 s1_list=list(0.405),
                 s2_list=list(0.305),
                 grid_scale_clever_covariate=0.01,
                 cvxr=T,
                 uniform_ratio=F)

  saveRDS(res, file = "071724_various/knT_071724_knownci_perc1_cvxrT_MLE_PenalUsual_smallPara_NoCollinear_0thOrder_withKnots&UpEnd_improve_2_" %+% n %+% "_" %+% seed %+% ".RDS")
}

#s1_list=list(0.205,0.205,0.205,0.405,0.405,0.405,0.605,0.605,0.605),
#s2_list=list(0.205,0.405,0.605,0.205,0.405,0.605,0.205,0.405,0.605),

# s1_list=list(0.205,0.405,0.405,0.605),
# s2_list=list(0.205,0.205,0.605,0.605),


# s1_list=list(0.205,0.205,0.405,0.405,0.505),
# s2_list=list(0.205,0.405,0.205,0.405,0.505),

# s1_list=list(0.405),
# s2_list=list(0.305),

#Small grid
# s1_list=list(0.105,0.105,0.305,0.305,0.205),
# s2_list=list(0.105,0.305,0.105,0.305,0.205),

#saveRDS(res, file = "041724_various/cn_041724_independent_usual_perc0.5_0_order_trueC_Trunc_cvxrT_ratioF_" %+% n %+% "_" %+% seed %+% ".RDS")

# seed_vals<-c(101,500,2000,4321,8688,10113,36899,45000,48001,50002)
# for (seed in seed_vals) {
#   n=100
#   res <- run_sim(n = n,
#                  seed = seed,
#                  type = "cn",
#                  nfolds=3,
#                  cv_choice="lambda",
#                  undersmooth_type="cv_all",
#                  basis_option="0_order",
#                  masking_option="no_masking",
#                  weight_option="no_tail_shrinkage",
#                  penalty_type="usual",
#                  undersmooth=F,
#                  censoring=F,
#                  perc_sample=1,
#                  perc_sample_basis=1,
#                  s1_list=list(0.2,0.4,0.4,0.6),
#                  s2_list=list(0.2,0.2,0.6,0.6),
#                  grid_scale_clever_covariate=0.01)
#
#   saveRDS(res, file = "030824_various/cn_031324_0.2corr_usual_full_0_order_rightpoint4_rightSubsampleAndPenalty_rightIntegration" %+% n %+% "_" %+% seed %+% ".RDS")
# }
#
#












#TODO: 1. tail shrinkage (put little emphasis on tail point) + flex_tail(tail will be fitted very flexibly). Both to ensure the
#           tail points have very limited effects for the middle but can still pool their information if strong. (DONE)
#     2. when subsampling, for the points greater than 0.7, force to add a points like 0.7 or use it to do the interval
#         dividing and then subsampling in each subinterval.(DONE)


#.   4. Try a different DGP other than a normal that centered around 0.5,0.5, and see what will be the bias curve looks like. I.E the observed whale shape of bias
#       may only be a special case to the normal DGP involved here. (DONE)
#.      Note in the file cl_022924 it is correlated dgp with 0.7corr. In addition, in both cl_022924, cl_030524, the propensity is calculated rather than estimated.

#.    5. use the estimated propensity scores.

#     3. Think about undersmoothing, maybe undersmoothing put too much emphasis on the tail region. If there is anyway for the centralic undersmoothing?
#         Try under smoothing with the masking OR shrinkage and/or with flex tail.
#     Q: To what degree the undersmoothing should be?
#         Flipping really help? more stable lasso's benefit? If flipping, how to get the update done


##Later time point, get rid of the 0 order indicators but get the first order splines to get extrapolation
#Old points are:
# s1_list=list(0.2,0.2,0.8,0.8),
# s2_list=list(0.2,0.8,0.2,0.8),

##flipping HAL fit


##Some regularization at the targeting step, more penalization, maybe more penalization on \epsilon_4 then \epsilon_1

#Targeting is for solving the EIC, rely on the aysmptotic behavior. But finite sample bias variance trade off may need a more
#regularized targeting, like slensky's method, we move the epsilon till the point where the change in variance is dominating the change in bias.

#Different versions:
#Weight 1 c(0.005,0.03,0.03,0.06)
#Weight 2 c(0.4,0.8,0.8,1)


#The point 0.2 0.8 for example, the hazards around there is only fitted very minimally fitted, but was used by nearly every point. So the targeting may hurt due to finite sample.
#We should only target for the point where the hazards is well estimated and hope they can bring their good effect to the edgy point like (0.2,0.8)
#Investigating:
#what info is needed in the estimation
#What info is needed in the targeting
#What info is really nailed in fitting

#May be just need to target center region
#New 4 points are:
# s1_list=list(0.2,0.4,0.4,0.6),
# s2_list=list(0.2,0.2,0.6,0.6),

#New 6 points are:
# s1_list=list(0.2,0.4,0.4,0.6,0.3,0.5),
# s2_list=list(0.2,0.2,0.6,0.6,0.5,0.3),

# x<-cbind(data_full$tt1,data_full$tt2)
# cl<-kmeans(x,4)
# plot(x, col = cl$cluster)
# points(cl$centers, col = 1:2, pch = 8, cex = 2)

# plot(data_full$t1,data_full$t2)
# plot(data_full$tt1,data_full$tt2)
