
make.dat <- function(setup.file.name, releases, recaptures, age.dependent=F,only.harvested=F,HM
                     ,eM=T,pM,sM,Mb # natural mortality settings
                     ,eF=T,pF,sF,Fb # fishing mortality settings
                     ,eT=T,pT,sT,Tb # tagging mortality settings
                     ,combine.Hrr_and_Rrr=F # combine Harvest & Release reporting rates
                     ,eHrr=F,pHrr,sHrr,Hrrb # Harvest reporting rates
                     ,eRrr=F,pRrr,sRrr,Rrrb # Release reporting rates
                     ,combine.HRS_and_RRS=F # combine Harvest & Release Retention Survival
                     ,eHRS=F,pHRS,sHRS,HRSb #  Harvest Retention Survival
                     ,eRRS=F,pRRS,sRRS,RRSb #  Release Retention Survival
                     ,incomplete.mix=F,pNonMix,sNonMix,NonMixb #  Non-Mixing
                     ,combine.HS_and_RS=T # combine Harvest & Release Selectivity
                     ,eHS=T,sHS,HSb # age-depdendent Harvest Selectivity
                     ,eRS=F,sRS,RSb # age-depdendent Release Selectivity
                     ,eMS=T,aMS,sMS # age-depdendent M Selectivity
                     ,use.L.const=T # Likelihood Constant
                     #,open.dat_file=T
){
  
  if(is.null(releases$y_rel) | is.null(recaptures$y_rec)) {
    releases$y_rel <- .as.year(releases$rel.date)
    recaptures$y_rel <- .as.year(recaptures$rel.date)
    recaptures$y_rec <- .as.year(recaptures$rec.date)
  }
  
  #' HM: hooking mortality per recapture year, only used if only.harvested is set FALSE
  #' eM: estimate natural Mortality
  #' pM: beginning years of different natural Mortality (M) periods (by default: one period for entire series)
  #' sM: starting value for natural mortality estimates (default value is 0.15)
  
  #   setup.file.name <- "test"
  #   age.dependent <- T
  #   only.harvested <- F
  if(is.null(releases$age) | is.null(recaptures$age)) {
    releases$age <- 1
    recaptures$age <- 1
    age.dependent <- F
    warning('releases and/or recaptures without age information! running age-independent model!')
  }
  
  if(!age.dependent) {
    releases$age <- 1
    recaptures$age <- 1
  }
  
  ages <- min(releases$age):max(recaptures$age)
  
  ########################### editing input data
  #### run release data:
  total_rel <- ddply(releases,c("y_rel","age"),function(x)c(n=nrow(x)))
  y_rel <- min(total_rel$y_rel):max(total_rel$y_rel)
  
  total_rel_list <- list()
  for(age in ages){
    m <- rep(0,length(y_rel))
    for(j in 1:length(y_rel)){
      i <- which(total_rel$y_rel == y_rel[j] & total_rel$age == age)
      if(length(i) > 0) m[j] <- total_rel$n[i]
      nm <- paste0("Total_Releases",.switch.if(age.dependent,paste0("_of_Age_Class_",which(ages==age)),""))
      total_rel_list[[nm]] <- m
    }
  }
  total_rel_list
  
  
  #### run recapture data:
  if(is.null(recaptures$harvested)){
    recaptures$harvested <- T
  }else{
    if(only.harvested){
      if(any(!recaptures$harvested)) warning('catch-release data treated as harvested!')
      recaptures$harvested <- T
    }
  }
  
  total_rec <- ddply(recaptures,c("y_rel","y_rec","age,harvested"),function(x)c(n=nrow(x)))
  y_rec <- min(total_rec$y_rec):max(total_rec$y_rec)
#   head(total_rec)
  
  ### create total_rec tables for harvested and catch & release sets
  total_rec_list <- list()
  for(harvested in c(T,F)){
    for(age in ages){  
      m <- matrix(nrow = length(y_rel),ncol=length(y_rec),-1)
      for(x in 1:length(y_rec)){
        for(y in 1:length(y_rel)){
          i <- which(total_rec$age == age & total_rec$y_rel == y_rel[y] & total_rec$y_rec == y_rec[x] & total_rec$harvested == harvested) 
          #         total_rec[i,]
          if(length(i) == 0 & y_rel[y] <= y_rec[x]) m[y,x] <- 0
          if(length(i) > 0) m[y,x] <- total_rec$n[i]
        }
        if(only.harvested & !harvested) m[,] <- -1
        nm <- paste0(.switch.if(harvested,"Harvest_","Release_"),"Tag_Recovery_Data",.switch.if(age.dependent,paste0("_for_Age_Class_",which(ages==age)),""))
        total_rec_list[[nm]] <- m
      }
    }
  }
  total_rec_list
  
  #### finished preediting data
  
  
  ################################## start writing input data file:
  
  fn0 <- paste0(setup.file.name,".dat")
  fn <- gsub('.dat.dat','.dat',fn0)
  zz <- file(fn, open = "wt")
  sink(zz)
  .catn('# IRATE Version 1.0')
  .catn(paste('# created by IRATER Version 1.0 -',date()))
  .catn(paste('# Model Type:',.switch.if(age.dependent,"Age-dependent","Age-Independent")))
  .catn(.truefalse2onekey(age.dependent))
  .catn('# Tag Recovery Type 0=har only, 1=har/rel')
  .catn(.switch.if(only.harvested,0,1))
  .catn('# Release Start & End Years') # tag release start and end years
  .catn(paste(range(y_rel),collapse=" "))
  .catn('# Recapture Start & End Years') # tag recapture start and end years
  .catn(paste(range(y_rec),collapse=" "))
  .catn('# Number of Ages') # number of age calsses
  .catn(length(ages))
  .catn('# First Age Class') # first age class
  
  .catn(min(total_rel$age))
  for(nm in names(total_rel_list)){
    .catn(paste('#', gsub('_'," ",nm))) # print total tag releases (per year and age class)
    m <- total_rel_list[[nm]]
    .catn(paste(m,collapse=" "))
  }
  
  for(nm in names(total_rec_list)){ # print both Harvest and Release Tag Recoveries
    .catn(paste('#', gsub('_'," ",nm))) # print total tag recaptures (per year and age class)
    m <- total_rec_list[[nm]]
    for(i in 1:nrow(m)) .catn(paste(m[i,],collapse=" "))
  }
  
  
  ################## Mortality settings:
  
  ### hooking mortality output: (only used in harvest+catch-and-release tagging studies)
  .catn('# Hooking Mortality')
  if(only.harvested){
    HM <- rep(-1,length(y_rec))
  }else{
    warning('setting hooking-mortality (HM) to default fixed value 0.09')
    if(missing(HM)) HM <- rep(0.09,length(y_rec))
  }
  .catn(paste(HM,collapse=" "))
  
  ### natural mortality output: 
  .catn('# Natural Mortality')
  .catn('# Estimate M  -1 = No, >=1 = Yes')
  .catn(.truefalse2onekey(eM))
  .catn('# Number of  M Periods')
  if(missing(pM)) pM <- min(y_rec) # default: assuming only one period of natural mortality
  .catn(paste(length(pM),collapse=" "))
  .catn('# M Periods Starting Years') # different periods of natural mortality?
  .catn(paste(pM,collapse=" "))
  .catn('# M Period Starting Values')
  if(missing(sM)) sM <- rep(0.15,length(pM)) # starting values of natural mortality
  .catn(paste(sM,collapse=" "))
  
  # set mixing settings:
  .catn('# Incomplete Mixing  -1 = No, >=1 = Yes')
  .catn(.truefalse2onekey(incomplete.mix))
  
  # fishing mortality settings:
  .catn('# Fishing Mortality')
  .catn('# Estimate F -1 = No, >=1 = Yes')
  .catn(.truefalse2onekey(eF))
  .catn('# Number of F Periods')
  if(missing(pF)) pF <- min(y_rec):max(y_rec) # default: assuming annual estimates of fishing mortality
  .catn(paste(length(pF),collapse=" "))
  .catn('# Place Holder Follows')  # either place holder and set to 0 or .catn('# Original') with a value
  .catn(paste(rep(0,2),collapse=" "))
  .catn('# F Period Starting Years')
  .catn(paste(pF,collapse=" "))
  .catn('# F Period Starting Values')
  if(missing(sF)) sF <- rep(0.1,length(pF)) # starting values of fishing mortality
  .catn(paste(sF,collapse=" "))
  
  # tagging mortality settings: (only for harvest and catch-release runs)
  .catn('# Tagging Mortality')
  .catn('# Estimate F-Tag  -1 = No, >=1 = Yes')
  if(only.harvested){
    if(!missing(eT)) warning('Release Reporting Rate-settings (eRrr & sRrr) disregarded as combine.Hrr_and_Rrr was set T')
    eT <- F
  }
  .catn(.truefalse2onekey(eT))
  .catn('# Number of F-Tag Periods')
  if(missing(pT)) pT <- min(y_rec)  # default: assuming single tagging mortality period
  .catn(paste(length(pT),collapse=" "))
  .catn('# Place Holder Follows')  # either place holder and set to 0 or .catn('# Original') with a value
  .catn(paste(rep(0,2),collapse=" "))
  .catn('# F-Tag Period Years')
  .catn(paste(pT,collapse=" "))
  .catn('# F-Tag Period Starting Values')
  if(missing(sT)) sT <- rep(0.1,length(pT)) # starting values for tagging mortality
  .catn(paste(sT,collapse=" "))
  
  
  
  ################## Reporting Rate settings:
  
  # Reporting Rate settings:
  .catn('# Combine L1 and L2') # equivalent to "Use Harvest Reporting Rate"; combine Harvest and Release Reporting Rate
  .catn(.truefalse2onekey(combine.Hrr_and_Rrr)) # always set -1, not found in JianADCR.dat; apparently a mistake in the file
  
  # Harvest Reporting Rate settings:
  .catn('# Harvest Reporting Rate') # L1
  .catn('# Estimate Harvest Reporting Rate -1 = No, >=1 = Yes')
  .catn(.truefalse2onekey(eHrr))
  .catn('# Number of Harvest Report Periods')
  if(missing(pHrr)) pHrr <- min(y_rec)  # default: assuming single Harvest reporting rate period
  .catn(paste(length(pHrr),collapse=" "))
  .catn('# Harvest Report Periods Years')
  .catn(paste(pHrr,collapse=" "))
  .catn('# Harvest Report Periods Starting Values')
  if(missing(sHrr)) sHrr <- rep(0.43,length(pHrr)) # starting values for Harvest reporting rates
  .catn(paste(sHrr,collapse=" "))
  
  # Release Reporting Rate settings:
  .catn('# Release Reporting Rate')
  .catn('# Estimate Release Reporting Rate -1 = No, >=1 = Yes')
  if(combine.Hrr_and_Rrr){
    if(!missing(eRrr)) warning('Release Reporting Rate-settings (eRrr & sRrr) disregarded as combine.Hrr_and_Rrr was set T')
    eRRS <- F
  }
  .catn(.truefalse2onekey(eRrr))
  .catn('# Number of Release Report Periods')
  if(missing(pRrr)) pRrr <- min(y_rec)  # default: assuming single Release reporting rate period
  .catn(paste(length(pRrr),collapse=" "))
  .catn('# Release Report Periods Years')
  .catn(paste(pRrr,collapse=" "))
  .catn('# Release Report Periods Starting Values')
  if(missing(sRrr)) sRrr <- rep(0.43,length(pRrr)) # starting values for Release reporting rates
  .catn(paste(sRrr,collapse=" "))
  
  
  
  ################## Retention-Survival settings:
  
  # Retention-Survival settings:
  .catn('# Combine phi1 and phi2') # equivalent to "Use Harvest RS"; combine Harvest and Release Retention Survival
  .catn(.truefalse2onekey(combine.HRS_and_RRS)) 
  
  # Harvest Retention Survival settings (phi):
  .catn('# Retention-Survival (phi)')
  .catn('# Estimate Harvest phi? -1= no; 1=yes')
  .catn(.truefalse2onekey(eHRS))
  .catn('# Number of phi-Harvest Periods')
  if(missing(pHRS)) pHRS <- min(y_rec)  # default: assuming single phi-Harvest period
  .catn(paste(length(pHRS),collapse=" "))
  .catn('# phi-Harvest Periods Years')
  .catn(paste(pHRS,collapse=" "))
  .catn('# phi-Harvest Starting Values')
  if(missing(sHRS)) sHRS <- rep(1,length(pHRS)) # starting values for phi-Harvest
  .catn(paste(sHRS,collapse=" "))
  
  # Release Retention Survival settings (phi):
  #   .catn('# Release Retention-Survival (phi)')
  .catn('# Estimate Release phi? -1= no; 1=yes')
  if(combine.HRS_and_RRS){
    if(!missing(eRRS)) warning('Release Retention-Survival-settings (eRRS & sRRS) disregarded as combine.HRS_and_RRS was set T')
    eRRS <- F
  }
  .catn(.truefalse2onekey(eRRS))
  .catn('# Number of phi-Release Periods')
  if(missing(pRRS)) pRRS <- min(y_rec)  # default: assuming single phi-Release period
  .catn(paste(length(pRRS),collapse=" "))
  .catn('# phi-Release Periods Years')
  .catn(paste(pRRS,collapse=" "))
  .catn('# phi-Release Starting Values')
  if(missing(sRRS)) sRRS <- rep(1,length(pRRS)) # starting values for phi-Release
  .catn(paste(sRRS,collapse=" "))
  
  ## Non-Mixing settings:
  .catn('# Number of Non-Mixing Periods')
  if(missing(pNonMix)) pNonMix <- min(y_rec)  # default: assuming single Non-Mix period
  .catn(paste(length(pNonMix),collapse=" "))
  .catn('# Non-Mixing Periods Years')
  .catn(paste(pNonMix,collapse=" "))
  .catn('# Non-Mixing Starting Values')
  if(missing(sNonMix)) sNonMix <- rep(0,length(pNonMix)) # starting values for phi-Release
  .catn(paste(sNonMix,collapse=" "))
  
  
  ################## Selectivity settings:
  
  # Selectivity settings (only applying in age-dependent models)
  .catn('# Age-dependent Selectivity')
  .catn('# Use harvest selectivities for both')
  #   .catn('Combine Harvest and Release selectivities -1 = No 1 = Yes') # Use harvest selectivities for both?
  .catn(.truefalse2onekey(combine.HS_and_RS)) 
  
  # Harvest Selectivity settings:
  .catn('# Estimate age-dependent Harvest Selectivity? -1= no; 1=yes') 
  .catn(.truefalse2onekey(eHS))
  .catn('# Harvest Selectivity Starting Values')
  if(missing(sHS)) sHS <- rep(0.5,length(ages)) # starting values for Harvest Selectivity
  .catn(paste(sHS,collapse=" "))
  
  # Release Selectivity settings:
  .catn('# Selectivity - Release')
  .catn('# Estimate age-dependent Release Selectivity? -1= no; 1=yes')
  if(combine.HS_and_RS){
    if(!missing(eRS)) warning('Release Selectivity-settings (eRS & sRS) disregarded as combine.HS_and_RS was set T')
    eRS <- F
  }
  .catn(.truefalse2onekey(eRS))
  .catn('# Release Selectivity Starting Values')
  if(missing(sRS)) sRS <- rep(0.5,length(ages)) # starting values for Release Selectivity
  .catn(paste(sRS,collapse=" "))
  
  # Natural Mortality Selectivity settings:
  .catn('# M selectivity')
  .catn('# Estimate M Selectivity? -1= no; 1=yes') # age-dependent Natural Mortality
  .catn(.truefalse2onekey(eMS))
  .catn('# Number of M Ages')
  if(missing(aMS)) aMS <- 1:length(ages)  # default: assuming unique Natural Mortality for each age class
  .catn(paste(length(aMS),collapse=" "))
  .catn('# M Beginning Ages')
  .catn(paste(aMS,collapse=" "))
  .catn('# M Selectivity Starting Values')
  if(missing(sMS)) sMS <- rep(1,length(aMS)) # starting values for age-dependent Natural Mortality
  .catn(paste(sMS,collapse=" "))
  
  # Likelihood settings
  .catn('# Use Likelihood Constant -1=no; 1=year')
  .catn(.truefalse2onekey(use.L.const))
  
  
  
  if(missing(Fb)) {
    Fb <- .switch.if(eF,c(-30,1.6),c(0,1))
  }else{
    if(!eF) warning('ignoring parameter boundary settings for Fishing Mortality (Fb), as F will not be estimated!')
  }
  
  if(missing(Mb)){
    Mb <- .switch.if(eM,c(-30,1.6),c(0,1))
  }else{
    if(!eM) warning('ignoring parameter boundary settings for Natural Mortality (Mb), as M will not be estimated!')
  }
  
  if(missing(Tb)){
    Tb <- .switch.if(eT,c(-30,1.6),c(0,1))
  }else{
    if(!eT) warning('ignoring parameter boundary settings for Hooking Mortality (Tb), as T will not be estimated!')
  }
  
  if(missing(Hrrb)){
    Hrrb <- .switch.if(eHrr,c(-30,1.6),c(0,1))
  }else{
    if(!eHrr) warning('ignoring parameter boundary settings for Harvest reporting rates (Hrrb), as Hrr will not be estimated!')
  }
  
  if(missing(Rrrb)){
    Rrrb <- .switch.if(eRrr,c(-30,1.6),c(0,1))
  }else{
    if(!eRrr) warning('ignoring parameter boundary settings for Release reporting rates (Rrrb), as Rrr will not be estimated!')
  }
  
  if(missing(HRSb)){
    HRSb <- .switch.if(eHRS,c(-30,1.6),c(0,1))
  }else{
    if(!eHRS) warning('ignoring parameter boundary settings for Harvest Retention-Survival (HRSb), as HRS will not be estimated!')
  }
  
  if(missing(RRSb)){
    RRSb <- .switch.if(eRRS,c(-30,1.6),c(0,1))
  }else{
    if(!eRRS) warning('ignoring parameter boundary settings for Release Retention-Survival (RRSb), as RRS will not be estimated!')
  }
  
  if(missing(HSb)){
    HSb <- .switch.if(eHS,c(-30,1.6),c(0,1))
  }else{
    if(!eHS) warning('ignoring parameter boundary settings for Harvest Selectivity (HSb), as HS will not be estimated!')
  }
  
  if(missing(RSb)){
    RSb <- .switch.if(eRS,c(-30,1.6),c(0,1))
  }else{
    if(!eRS) warning('ignoring parameter boundary settings for Release Selectivity (RSb), as RS will not be estimated!')
  }
  
  if(missing(NonMixb)){
    NonMixb <- .switch.if(incomplete.mix,c(-30,1.6),c(0,1))
  }else{
    if(!incomplete.mix) warning('ignoring parameter boundary settings for Incomplete Mixing (NonMixb), as no incomplete mixing was assumed!')
  }
  
  # Parameter Boundary settings (Controls):
  .catn('# Lower and Upper Boundary Values') # only considered if parameter is set to be estimated
  .catn(paste(Fb,collapse=' '))
  .catn(paste(Mb,collapse=' '))
  .catn(paste(Tb,collapse=' '))
  .catn(paste(Hrrb,collapse=' '))
  .catn(paste(Rrrb,collapse=' '))
  .catn(paste(HRSb,collapse=' '))
  .catn(paste(RRSb,collapse=' '))
  .catn(paste(HSb,collapse=' '))
  .catn(paste(RSb,collapse=' '))
  .catn(paste(NonMixb,collapse=' '))
  
  sink()
  #if(open.dat_file) system(paste('gedit',fn))
}



.catn <- function(x,...){
  cat(paste0(x,"\n"),...)
}

.truefalse2onekey <- function(x){
  if(x){1}else{-1}
}

.as.year <- function(x){
  as.numeric(format(as.Date(x),"%Y"))
}

.switch.if <- function(x,a,b){
  if(x){
    y <- a
  }else{
    y <- b
  }
  return(y)
}

