# prepare data and run missing data model
# constrain on 3, 4, 11

rm(list=ls())

wkdir="C:/Users/JZhang/Dropbox/Missing_Data"   # home
wkdir="C:/Users/jz250/Dropbox/Missing_Data"   # duke

setwd(wkdir)

load("../PPMI/Clean_Data/data0123.RData")

#####################################
# four file name  MDS_UPDRS_surv, event_table, screen_post, ID_surv

head(UPDRS_valid)

UPDRS_temp=UPDRS_valid
########################################
# delete several column

UPDRS_temp$INFODT=NULL
head(UPDRS_temp)
dim(UPDRS_temp)  #3309   68

dim(UPDRS_temp[,c(3:61)])  # 4921 59
##########################
# force level start from 1
UPDRS=UPDRS_temp[,c(3:61)]+1
head(UPDRS)
dim(UPDRS)
#######################
# test NA
test_NA=rep(0, nrow(UPDRS))
sum(UPDRS[1,])
for(i in 1: nrow(UPDRS)){
  if(is.na(sum(UPDRS[i,]))) test_NA[i]=1
  
}

test_NA-UPDRS_temp$rr
sum(abs(test_NA-UPDRS_temp$rr))   # 0 good
#UPDRS_temp[which(test_NA-UPDRS_temp$rr !=0 ),]  # patno 4070 has baseline NA
# UPDRS_temp[UPDRS_temp$PATNO==4070,]

rr_obs=UPDRS_valid$rr


################################
# replace NA as -1, otherwise STAN will not recognize


for(i in 1: nrow(UPDRS)){
  if(is.na(sum(UPDRS[i,]))) UPDRS[i,]= -1
  
}

head(UPDRS)
dim(UPDRS)   # 3309   59
Y_ordi_part1=UPDRS[,c(1:13)]
Y_ordi_part2=UPDRS[,c(14:26)]
Y_ordi_part3=UPDRS[,c(27:59)]

head(Y_ordi_part1)
head(Y_ordi_part2)
head(Y_ordi_part3)
#######################
#use 
p1=length(which(Y_ordi_part1[,3]==1))/nrow(UPDRS); p1   # 0.684
a1=log(p1/(1-p1)); a1  # 0.771718

p2=length(which(Y_ordi_part2[,4]==1))/nrow(UPDRS);p2  # 0.548

a2=log(p2/(1-p2)); a2   # 0.19

p3=length(which(Y_ordi_part3[,11]==1))/nrow(UPDRS);p3  # 0.3348

a3=log(p3/(1-p3)); a3   # -0.68

# a0=c(0.7,0.15, -0.65)

 a0=c(1.5,0.8, -0.5)




sum(UPDRS_valid$rr)   # 158, 158 missing


summary(unlist(UPDRS))   # no na good

num_subject=length(unique(UPDRS_valid$PATNO)); num_subject   # 415

num_obs= nrow(UPDRS_valid); num_obs   # 3309

num_part1=ncol(Y_ordi_part1); num_part1   # 13

num_part2=ncol(Y_ordi_part2); num_part2   # 13

num_part3=ncol(Y_ordi_part3); num_part3   # 33



num_ordi=5

subj_long=UPDRS_valid$subject
length(subj_long)   # 3309

unique(subj_long)-seq(1:415)   # all 0 good




age_norm=(screen_post$age- mean(screen_post$age))/sd(screen_post$age)
head(screen_post)
gender_subj=screen_post$gender

#######################################
# tee and time use year
head(event_table_valid)
tee=event_table_valid$tee

head(tee)  # make sure in year

head(UPDRS_valid[,c(1:3,65:69)],16)

time_obs=UPDRS_valid$time/12

summary(time_obs)   # 0.0000  0.3333  1.0000  1.6800  3.0000  6.4167 
head(time_obs,15)
tail(tee)
tail(time_obs) # good
hist(tee)
summary(tee)
table(event_table_valid$status)


# 178 237 

summary(time_obs)
#  0.0000  0.3333  1.0000  1.6800  3.0000  6.4167 
summary(tee)
# 0.1667  0.5000  2.1667  2.6994  4.9167  6.4167

###############################################
# piece-wise hazard setting
tau=seq(0, 7, by=1); tau  # 5 cut value use ceiling(max(tee)) as end point
num_tau<-length(tau) ; num_tau   # 5
n.intv= length(tau)-1 ; n.intv    # 4 interval


# D_tau is the tau difference
# D_tau <- tau[2:n.cutpt] - tau[1:n.intl] # use this vector to assure the interval may not equal, otherwise use 1.5 is enough
D_t<- matrix(NA, num_subject, n.intv)
seq_tau<-matrix(0, num_subject, num_tau)
tee_id<- rep(0, num_subject)
Ind <-matrix(NA, num_subject, n.intv)    # indicate tee in which interval of tau
#teeV_id<-matrix(0, num_subj, n.intl) # use cum_ind to inicate if interval within tee
for(i in 1: num_subject){
  for(k in 1: n.intv){
    Ind[i,k]=ifelse(tee[i]>=tau[k] & tau[k+1]>tee[i], 1, 0)  # tauR1 start from 0
    D_t[i,k]<- ifelse(tee[i]>tau[k], ifelse(tee[i]>=tau[k] & tau[k+1]>tee[i], tee[i]-tau[k], tau[k+1]-tau[k]),0) 
    seq_tau[i,k+1]<- ifelse(tee[i]>tau[k], ifelse(tee[i]>=tau[k] & tau[k+1]>tee[i], tee[i], tau[k+1]),tee[i]) 
    
  }
  if(sum(Ind[i,])==0) tee_id[i]=n.intv else tee_id[i]=which(Ind[i,]==1)
}


# verify
head(Ind)
head(tee)
head(tee_id) #good
head(seq_tau)   # good

tail(tee)
tail(tee_id)
tail(seq_tau)


num_pw=n.intv; num_pw   # 7



library(rstan)


data <- list(num_obs=num_obs, num_subject=num_subject, num_part1=num_part1,num_part2=num_part2,num_part3=num_part3,
              num_ordi=num_ordi, subj_long=subj_long,  Y_ordi_part1=Y_ordi_part1, Y_ordi_part2=Y_ordi_part2,
             Y_ordi_part3=Y_ordi_part3,
             a0=a0,  time_obs=time_obs,  age_norm=age_norm,  gender_subj=gender_subj,
             tee=tee, event=event_table_valid$status,  rr_obs=rr_obs, num_pw=num_pw, subj_pw_ind=tee_id, tee_pw=seq_tau)


pars <- c("beta0","beta1", "alpha", "a_ordi_part1", "a_ordi_part2", "a_ordi_part3",
           "b_ordi_part1", "b_ordi_part2", "b_ordi_part3", 
               "Omega", "Var_U",  "Var_e",  "w", "eta",
               "sd_U", "sd_e",
              "gam",  "nu","h0" )


pars_ConsE <- c("beta1", "alpha", "a_ordi_part1", "a_ordi_part2", "a_ordi_part3",
          "b_ordi_part1", "b_ordi_part2", "b_ordi_part3", 
          "Omega", "Var_U",   "w", "eta",
          "sd_U",       "gam",  "nu","h0" )

inits01 <- list( U=matrix(0.1, num_subject, 6),  Omega= diag(6), Var_e=rep(1,3),
                     Var_U=rep(1, 6),   ee=matrix(0, 3, num_obs), 
                     beta0=rep(0,3),  beta1=rep(0, 3), alpha=rep(0.1,3), 
                     a_random1=rep(0.9, num_part1),  a_random2=rep(0.8, num_part2), a_random3=rep(1, num_part3),
                     b_random1= rep(0.3, num_part1), b_random2= rep(0.4, num_part2), b_random3= rep(0.2, num_part3),
                     delta1=matrix(1, nrow=num_part1, ncol=num_ordi-2),  delta2=matrix(1, nrow=num_part2, ncol=num_ordi-2), 
                     delta3=matrix(1, nrow=num_part3, ncol=num_ordi-2), 
                     gam=-0.1, nu=rep(0.3,3), h0=rep(0.005, num_pw), w=-8, eta=rep(1,3) )

inits01 <- list(c1=inits01)



inits02 <- list( U=matrix(0.3, num_subject, 6),  Omega= diag(6), Var_e=rep(2,3),
                 Var_U=rep(2, 6),   ee=matrix(-0.1, 3, num_obs), 
                 beta0=rep(-1,3),  beta1=rep(0.5, 3), alpha=rep(0.2,3), 
                 a_random1=rep(1, num_part1),  a_random2=rep(1, num_part2), a_random3=rep(0.8, num_part3),
                 b_random1= rep(0.1, num_part1), b_random2= rep(0.2, num_part2), b_random3= rep(0.3, num_part3),
                 delta1=matrix(0.5, nrow=num_part1, ncol=num_ordi-2),  delta2=matrix(1.1, nrow=num_part2, ncol=num_ordi-2), 
                 delta3=matrix(1.5, nrow=num_part3, ncol=num_ordi-2), 
                 gam=0.1, nu=rep(-0.2,3), h0=rep(0.001, num_pw), w=-2, eta=rep(0.5,3) )


inits02 <- list(c1=inits02)


inits01_ConsE <- list( U=matrix(0.1, num_subject, 6),  Omega= diag(6), 
                 Var_U=rep(1, 6),   ee=matrix(0, 3, num_obs), 
                   beta1=rep(0, 3), alpha=rep(0.1,3), 
                 a_random1=rep(0.9, num_part1),  a_random2=rep(0.8, num_part2), a_random3=rep(1, num_part3),
                 b_random1= rep(0.3, num_part1), b_random2= rep(0.4, num_part2), b_random3= rep(0.2, num_part3),
                 delta1=matrix(1, nrow=num_part1, ncol=num_ordi-2),  delta2=matrix(1, nrow=num_part2, ncol=num_ordi-2), 
                 delta3=matrix(1, nrow=num_part3, ncol=num_ordi-2), 
                 gam=-0.1, nu=rep(0.3,3), h0=rep(0.005, num_pw), w=-8, eta=rep(1,3) )

inits01_ConsE <- list(c1=inits01_ConsE)


#############################################
model_file<-"./missing_11_real_3LVs.stan"
model_file_ConsE= "./missing_11_real_3LVs_ConsE.stan"

time0<-Sys.time()
fit1<- stan(file=model_file, data=data, pars=pars, init=inits01,  thin=1, chains=1, iter=3000, warmup=2000, seed=1234)
Sys.time()-time0 # 6.527301 hours

time0<-Sys.time()
fit1_ConsE<- stan(file=model_file_ConsE, data=data, pars=pars_ConsE, init=inits01_ConsE,  thin=1, chains=1, iter=3500, warmup=2000, seed=1234)
Sys.time()-time0 # 6.527301 hours


time0<-Sys.time()
fit2<- stan(file=model_file, data=data, pars=pars, init=inits02,  thin=1, chains=1, iter=3500, warmup=2500, seed=1234)
Sys.time()-time0 # 6.618393 hour


save(fit2, file = 'U:/missing_data/fit2')

pars_est= c("beta0", "beta1", "alpha",  "Omega", "sd_U",   "w", "eta",
            "sd_U", "sd_e",
            "gam",  "nu","h0" )

main_est=summary(fit2, pars=pars_est, probs=c(0.025,0.975))$summary


# main_est=summary(fit1_ConsE, pars=pars_est, probs=c(0.025,0.975))$summary

print(main_est[,c(1,3:5,7)],digits=3)

library(xtable)

xtable(main_est[,c(1,3:5)], digits=3)
pars_ordinal1=c("a_ordi_part1", "b_ordi_part1")
pars_a1=c("a_ordi_part1")
a1_est=summary(fit2, pars=pars_a1, probs=c(0.025,0.975))$summary
ordinal1_est=
xtable(ordinal1_est[,c(1,3:5)], digits=3)

a1_ordinal=matrix(as.numeric(a1_est[,1]),13, 4, byrow=T) 

xtable(a1_ordinal, digits=3)

b1