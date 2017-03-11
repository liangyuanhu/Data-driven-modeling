#----------------------------------------------------------------#
#                      sample code                               #
#                ENAR 2017 Spring Meeting                        #
# T3: Data-driven modeling techniques in medical decision making #
#----------------------------------------------------------------#

rm(list=ls())
library('reshape')
library('reshape2')
library('dplyr')
library('survival') 
library('zoo')
library('lattice')
library(iterators)
library(itertools)
library(foreach)   

#because of sensitivity of our data, we can provide a small subset of the data for illustration
datsamp<- read.csv("datsamp.csv")
#----------------------------------------#
# to visualize CD4% trajectory over time #
#----------------------------------------#
xyplot(cd4per ~ fptime/30 |factor(patient_id), data = datsamp, groups = as.factor(onarv),
       auto.key = list(text=c("not initiated","initiated"), cex = 2, columns=2),
       par.settings = list(superpose.symbol = list(col = 1, pch = c("o", "+"), cex = 2.5, lwd =2)), 
       xlab = list(label ="Follow-up time (in months)", cex = 2),
       aspect=0.5, layout = c(2,4), ylab= list(label ="cd4 percentage",cex=2 ), 
       scales=list(cex=2),
       subset = patient_id %in% c(3802,3818,4240,4251,4330, 4722, 5922, 25179))

#------------------------------------------------#
#  bootstrap patients with longitudinal records  #
#------------------------------------------------#
clusters  =  as.numeric ( as.character( factor(names(table(datsamp$patient_id)) ) ))
cd4set  =  c(99, seq(25,5,by = -5))
set.seed (2950)
index     =    sample(1:length(clusters), length(clusters), replace=TRUE)

aa        =    clusters[index]  #sampled id 

bb        =    table(aa)   # see how many times each id was sampled 

boot.samp   =  NULL

for(j in 1:max(bb) ){
  
  #if an id was sampled j times,  then select, j times, observations corresponding to that id 
  
  cc  =  data[data$patient_id %in% names(bb[bb %in% j]),]  
  
  for(k in 1:j){                                 
    
    boot.samp  =  rbind(boot.samp, cc)
    
  }
  
}

boot.samp = boot.samp[order(boot.samp$patient_id),]

dd1  =  table(boot.samp$patient_id) #no. of rows for each patient  
dd2  =  dd1/bb
dd3  =  rep(dd2, bb)

#use subjid, treat patients that are sampled x times as x subjects 
boot.samp$subjid  = rep(1:length(clusters), dd3)

#---------------------------------------------------------#
# function to create expanded regimen data
#---------------------------------------------------------#

ifun = function (x) {
  
  idata = filter(data, subjid == x)
  
  n = dim(idata)[1]
  
  fptime     =       idata$fptime   
  
  onarv_t    =       idata$onarv 
  
  cd4per     =       idata$cd4per
  
  last       =       length(cd4per)
  
  # cd4per_af_1 =      cd4per[2:last] #after 1st clinical visit
  
  mincd4      =      idata$mincd4 
  
  cd4permr     =     idata$cd4permr 
  
  
  #=======================================================================================================================
  #create multiple regimes 
  
  regime     = list()  #regimes followed at each follow up visit
  indicator  = list()  #indicators for each regime  (e.g., if R25 = 1, following d = 25 )
  
  
  if (n > 3) { onarv_lead_3  =  lead(onarv_t, 3);  onarv_lead_3[(n-2):n] = onarv_t[(n-2):n];
  
  onarv_lead_2  =  lead(onarv_t, 2);  onarv_lead_2[(n-1):n] = onarv_t[(n-1):n];
  
  onarv_lead_1  =  lead(onarv_t, 1);  onarv_lead_1[n] = onarv_t[n-1];
  
  onarv         =   onarv_t;
  
  #onarv[1:(n-3)]   =  onarv_lead_3[1:(n-3)] ;
  
  for (i in 1:n) {
    
    if (is.na (cd4per[i]) == 0 & sum(onarv_lead_1[i],onarv_lead_2[i], onarv_lead_3[i]) > 0 ) { onarv[i] = 1 }
    
  }
  
  onarv_lag_1     = lag(onarv, 1);
  
  }
  
  if (n < 4) { onarv  =   onarv_t; onarv_lag_1 = lag(onarv, 1) }
  
  
  if (i > 1) {
    
    # no missing baseline cd4, so no missing  "mincd4"
    # if (onarv_lag_1[i] == 0 & onarv[i] == 0 & is.na(mincd4[i]) == 1)  (regime[[i]] = cd4set[-1]) ;
    
    if (onarv_lag_1[i] == 0 & onarv[i] == 0 )      (regime[[i]] = cd4set[cd4set < mincd4[i]]);
    
    if (onarv_lag_1[i] == 0 & onarv[i] == 1 )      (regime[[i]] = regime[[i-1]] [ regime[[i-1]]>= cd4permr[i] ] ); #follow previous regimes where regimes >=current cd4%                                                          
    
    #once initiated, always on, always follow previous regimes 
    
    if (onarv_lag_1[i] == 1)                      (regime [[i]] = regime [[i-1 ]]  ); 
    
    indicator[[i]] = as.numeric(cd4set %in% regime[[i]]) 
    
  }
  
  
  
  regimes              =    as.data.frame( t(sapply(regime, '[', 1:6) ) )
  
  #maximum follow up length, if < n, then artifically censored 
  fl.max.len           =    max(apply(regimes,2,function(x) sum(!is.na(x) == 1 )))
  
  colnames(regimes)    =    c("reg1", "reg2", "reg3", "reg4", "reg5", "reg6") 
  
  indicators           =    as.data.frame( t(sapply(indicator, '[', 1:6) ) )
  
  colnames(indicators) =    c(99, 25, 20, 15, 10, 5)  
  
  idata_bind      =  cbind(select(idata, subjid, fptime),  indicators)
  
  idata_bind$fl.max.len   =  fl.max.len 
  
  #use reshape package function "melt" to get "molten" data at each follow up time 
  
  idatam    =    melt(idata_bind,id  =  c("subjid","fptime","fl.max.len"))
  
  idatam$art.C  =    ifelse(idatam$fl.max.len == n, 0, 1) 
  #convert characters to numeric --- as.character first, then as.numeric to get proper values
  #get regimes (numeric) 
  idatam$d = as.numeric(as.character(factor(idatam$variable)))   
  
  idata_1  =  filter(idatam, value == 1) 
  
  
  #define replicates; unique d and repeat times 
  if (nrow(idata_1) == 0) (dat_out = NULL)
  
  if (nrow(idata_1) > 0) 
    
  {
    
    reg_no  =  length(unique(idata_1$d));  #regimes followed 
    
    idata_1 = group_by(idata_1, d) %>% mutate(., n = n()); #each regime was followed n times 
    
    rep_times = unique(cbind(idata_1$d,idata_1$n));
    
    idata_1$replicate = rep(1:reg_no, rep_times[,2]);
    
    idata_1             =    select(idata_1, subjid, fptime, replicate, d, art.C);
    
    dat_out             =    left_join(idata_1, idata, by = c("subjid","fptime")) 
  }
  
  return(dat_out)
  
} 

IDs    =  unique(data1$subjid)
dat_expd =  foreach (x = IDs, .combine = 'rbind')%do% ifun(x) 
grpdat_expd  =  group_by(dat_expd, subjid)



#-----------------------------------------#
#   Use random forests for weight model   #
#-----------------------------------------#
#datsamp2 is the subset of datsamp with ART first time initiated

wmod_den <- tuneRF(x=datsamp2[,c("wazmr", "hazmr", "cd4mr", "male", "age", "education", "urban", "WHO", "CDC", "ethnic", "orphan", "periods")], 
                   y=datsamp2[,c("onarv")], mtryStart=1, ntreeTry=500, stepFactor=2, doBest=T)

wmod_fit = predict(wmod_den)

wmod_fit[ wmod_fit > 1] <- 1  ## correct for possible overflow by Random Forest due to machine accuracy
wmod_fit[ datsamp$onarv == 0] <- 1 - wmod_fit[datsamp$onarv == 0]

ids <- unique(datsamp2$patient_id)

ipcw_den <- NULL
for (i in ids)
{
  tmp <- datsamp2[datsamp2$patient_id==i,]
  mod_tmp <- wmod_fit[datsamp2$patient_id==i]
  ipcw_d <- mod_tmp
  tmp_d <- cumprod(ipcw_d)
  ipcw_den <- c(ipcw_den, tmp_d)
}
datsamp2$wt <- 1/ipcw_den  # cal denominator first, time varying 


#=========================================================================================================#
# impute missing outcome  
# imputation model:  probit(cd4per_12/100) = c0 + c1*cdcclass_0 + c2*cd4per_b + c3*onarv + c4*cdcclass_12
#=========================================================================================================#

gregmdat   =  group_by(regm_dat, by = subjid)

covX_g = gregmdat %>% summarise_each ( funs(tail(.,1)), cdcclass0, cd4per_b, onarv,cdcclass, cd4per_12, waz0, haz0, waz_12, haz_12)
size  = dim(covX_g)[1]

#convert cdcclass to numeric 
cdcclass0  = ifelse(covX_g[,2] == "A",1,ifelse(covX_g[,2] == "B", 2, ifelse(covX_g[,2] == "C",3, ifelse(covX_g[,2] == "N", 4, 5))))
cdcclass_12  = ifelse(covX_g[,5] == "A",1,ifelse(covX_g[,5] == "B", 2, ifelse(covX_g[,5] == "C",3, ifelse(covX_g[,5] == "N", 4, 5))))


cd4_12  =  as.numeric(unlist(covX_g[,6]))/100 
cd4_12[which(cd4_12 == 0)] = .01 

#inverse cummulative normal - inverse probit 
cd4_12[is.na(cd4_12) == 0] = qnorm(cd4_12[is.na(cd4_12) == 0] )

#missing indicator for cd412
R_cd4  =  ifelse(is.na(cd4_12) ==1, 0, 1 )  

#===============================================#
#z scores
waz_12  = as.numeric(unlist(covX_g[,9]))
haz_12  = as.numeric(unlist(covX_g[,10]))
R_waz   = ifelse(is.na(waz_12) ==1, 0, 1 )  
R_haz   = ifelse(is.na(haz_12) ==1, 0, 1 ) 

imp.dat = data.frame (cdc0 = cdcclass0, cd4per_b = as.numeric(unlist(covX_g[,3])), onarv = as.numeric(unlist(covX_g[,4])), 
                      cdc12= cdcclass_12, Y = cd4_12, R_cd4 = R_cd4)

imp.dat_waz = data.frame (cdc0 = cdcclass0, cd4per_b = as.numeric(unlist(covX_g[,3])), waz0 = as.numeric(unlist(covX_g[,7])),
                          onarv = as.numeric(unlist(covX_g[,4])), cdc12= cdcclass_12, Y = waz_12, R_waz = R_waz)

imp.dat_haz = data.frame (cdc0 = cdcclass0, cd4per_b = as.numeric(unlist(covX_g[,3])), haz0 = as.numeric(unlist(covX_g[,8])),
                          onarv = as.numeric(unlist(covX_g[,4])), cdc12= cdcclass_12, Y = haz_12, R_haz = R_haz)

X         =  as.matrix(cbind ( rep(1, size), imp.dat[,1:4]))
X_waz     =  as.matrix(cbind ( rep(1, size), imp.dat_waz[,1:5]))
X_haz     =  as.matrix(cbind ( rep(1, size), imp.dat_haz[,1:5]))


names(imp.dat) = c("cdc0", "cd4per_b", "onarv", "cdc12", "Y", "R_cd4")
names(imp.dat_waz) = c("cdc0", "cd4per_b", "waz0", "onarv", "cdc12", "Y", "R_waz")
names(imp.dat_haz) = c("cdc0", "cd4per_b",  "haz0"  ,"onarv", "cdc12", "Y", "R_haz")

fit.imp  =  lm (Y ~ cdc0 + cd4per_b + onarv + cdc12 , data = imp.dat, subset = (R_cd4==1) )

beta.cov  = vcov(fit.imp)

#imputation model residual 
tau      =  summary(fit.imp)$sigma


fit.imp_waz  =  lm (Y ~ cdc0 + cd4per_b + waz0 + onarv + cdc12 , data = imp.dat_waz, subset = (R_waz==1) )

fit.imp_haz  =  lm (Y ~ cdc0 + cd4per_b + haz0 + onarv + cdc12 , data = imp.dat_haz, subset = (R_haz==1) )


beta.cov_waz  = vcov(fit.imp_waz)

tau_waz      =  summary(fit.imp_waz)$sigma


beta.cov_haz  = vcov(fit.imp_haz)

tau_haz      =  summary(fit.imp_haz)$sigma

#------------------------------------------------------------------------#
#multiple imputation 200 times
#draw betas 

beta.samp  =  mvrnorm(1,mu = fit.imp$coef,Sigma = beta.cov)

Y.temp  =   X %*% beta.samp

Y.imp = pnorm(rnorm(size, Y.temp, tau))*100

#replace non-missing with observed values

Y.imp[imp.dat$R_cd4==1]  =  pnorm(imp.dat$Y[imp.dat$R_cd4==1])*100
sid =  unique(gregmdat$subjid)

imp.out  = data.frame (subjid = sid, cd4per_12 = Y.imp)


