# prepare data and run missing data model

rm(list=ls())

wkdir="C:/Users/JZhang/Dropbox/Missing_Data"   # home
wkdir="C:/Users/jz250/Dropbox/Missing_Data"   # duke

setwd(wkdir)

load("../PPMI/Clean_Data/data0118.RData")

#####################################
# four file name  MDS_UPDRS_surv, event_table, screen_post, ID_surv

head(MDS_UPDRS_surv)

UPDRS_temp=MDS_UPDRS_surv
########################################
# delete several column

UPDRS_temp$INFODT=NULL
head(UPDRS_temp)
dim(UPDRS_temp)  #4921   69

dim(UPDRS_temp[,c(3:61)])  # 4921 59
##########################
# force level start from 1
UPDRS=UPDRS_temp[,c(3:61)]+1


test_NA=rep(0, nrow(UPDRS))
sum(UPDRS[1,])
for(i in 1: nrow(UPDRS)){
  if(is.na(sum(UPDRS[i,]))) test_NA[i]=1
  
}
test_NA-UPDRS_temp$rr
sum(abs(test_NA-UPDRS_temp$rr))   # 1
UPDRS_temp[which(test_NA-UPDRS_temp$rr !=0 ),]  # patno 4070 has baseline NA
UPDRS_temp[UPDRS_temp$PATNO==4070,]
 ###############################################################
# need delete 4070's baseline
# UPDRS_valid
dim(MDS_UPDRS_surv)   # 4921 70
UPDRS_valid=MDS_UPDRS_surv[which(abs(test_NA-MDS_UPDRS_surv$rr) ==0 ),]

dim(UPDRS_valid)  # 4920 70
UPDRS_valid$INFODT=NULL


dim(UPDRS_valid)   # 4920   69

UPDRS=UPDRS_valid[,c(3:61)]+1

dim(UPDRS)  # 4920   59

head(UPDRS)


test_NA=rep(0, nrow(UPDRS))
# sum(UPDRS[1,])
for(i in 1: nrow(UPDRS)){
  if(is.na(sum(UPDRS[i,]))) test_NA[i]=1
  
}
sum(abs(test_NA-UPDRS_valid$rr) )   # 0 good

sum(UPDRS_valid$rr)   # 351, 351 missing

################################
# replace NA as -1, otherwise STAN will not recognize

for(i in 1: nrow(UPDRS)){
  if(is.na(sum(UPDRS[i,]))) UPDRS[i,]= -1
  
}

summary(unlist(UPDRS))   # no na good

num_subject=length(unique(UPDRS_valid$PATNO)); num_subject   # 421

num_obs= nrow(UPDRS_valid); num_obs   # 4920

num_item=ncol(UPDRS); num_item   # 59

num_ordi=5
subj_long=UPDRS_valid$subject
length(subj_long)
Y_ordi=UPDRS

length(which(Y_ordi[,1]==1)); length(which(Y_ordi[,1]==2)); length(which(Y_ordi[,1]==3));length(which(Y_ordi[,1]==-1))

p0=length(which(Y_ordi[,1]==1))/num_obs; p0   # 0.6132114

log(p0/(1-p0));   # 0.4608314

a0=0.5

age_subj=screen_post$age
head(screen_post)
gender_subj=screen_post$gender

#######################################
# tee and time use year
tee=event_table$tee/12
head(tee)
time_obs=UPDRS_valid$time/12
head(time_obs, 20)
event=event_table$status
head(event)


rr_obs=UPDRS_valid$rr
head(rr_obs,20)
tail(rr_obs, 20)
first_last_obs=UPDRS_valid$first_last
head(first_last_obs, 20)
tail(first_last_obs, 20)

library(rstan)


data <- list(num_obs=num_obs, num_subject=num_subject, 
             num_item=num_item, num_ordi=num_ordi, subj_long=subj_long,  Y_ordi=Y_ordi, 
             a0=a0,  time_obs=time_obs,  age_norm=age_subj,  gender_subj=gender_subj,
             tee=tee, event=event,  rr_obs=rr_obs, first_last_obs=first_last_obs)


pars <- c("beta", "alpha", "a_ordi", "b_ordi", 
               "Omega", "Var_U",  "Var_e",  "w", "eta",
               "sd_U", "sd_e",
              "gam",  "nu","h0" )

inits01 <- list(U=matrix(0.1, num_subject, 2),  Omega= diag(2), Var_e=1,
                     Var_U=rep(2, 2),   ee=rep(0, num_obs), 
                     beta=rep(0,2), alpha=0.1, 
                     a_random=rep(0.9, num_item), b_random= rep(0.3, num_item),
                     delta=matrix(1, nrow=num_item, ncol=num_ordi-2),   gam=-0.1, nu=0.3, h0=0.005, w=-8, eta=1 )

inits01 <- list(c1=inits01)



inits02 <- list(U=matrix(0.2, num_subject, 2),  Omega= diag(2), Var_e=2,
                     Var_U=rep(1, 2),  ee=rep(0.1, num_obs), 
                     beta=rep(0.5,2), alpha=0.2, 
                     a_random=rep(1, num_item), b_random= rep(0.5, num_item),
                     delta=matrix(0.5, nrow=num_item, ncol=num_ordi-2),  gam=0.1,  nu=0.2, h0=0.002, w=-5, eta=2 )

inits02 <- list(c1=inits02)


#############################################
model_file<-"./missing_11_real.stan"


time0<-Sys.time()
fit1<- stan(file=model_file, data=data, pars=pars, init=inits01,  thin=1, chains=1, iter=3500, warmup=2500, seed=1234)
Sys.time()-time0 # 6.527301 hours


print(fit1, digits=3)


time0<-Sys.time()
fit2<- stan(file=model_file, data=data, pars=pars, init=inits02,  thin=1, chains=1, iter=3500, warmup=2500, seed=1234)
Sys.time()-time0 # 2.031872 days


print(fit2, digits=3)

pars_est= c("beta", "alpha",  "Omega", "Var_U",  "Var_e",  "w", "eta",
            "sd_U", "sd_e",
            "gam",  "nu","h0" )

main_rst=summary(fit2, pars=pars_est, probs=c(0.025,0.975))$summary

library(xtable)
xtable(main_rst[,c(1,3:5)], digits=3)
