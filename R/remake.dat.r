remake.dat <- function(old.setup.name, new.setup.name="test", age.dependent=T,only.harvested = T,...){
  
  k <- read.dat(old.setup.name) # read example data set (with own function)
  df <- c()
  y_rel <- k$Release_Start_End_Years[1]:k$Release_Start_End_Years[2]
  y_rec <- k$Recapture_Start_End_Years[1]:k$Recapture_Start_End_Years[2]
  
  
  if(any(grep('Age_Class_',names(k)))){   ### check if age dependent:
    for(i in grep('Age_Class_',names(k)))
    {
      print(names(k)[i])
      a0 <- k[[names(k)[i]]]
      a0[a0 == -1] <- NA
      if(grepl('Total_Releases',names(k)[i])){
        add <- data.frame(n=a0$n,y_rel=a0$year,y_rec=NA)
        add$file <- names(k)[i]
        add$age <- as.numeric(substr(add$file,nchar(add$file),nchar(add$file)))
        add$file <- substr(add$file,1,(nchar(add$file)-1))
        df <- rbind(df,add)
      }else{
        for(y in 1:nrow(a0)){
          for(x in 1:ncol(a0)){
            add <- data.frame(n=a0[y,x],y_rel=y_rel[y],y_rec=y_rec[x])
            add$file <- names(k)[i]
            add$age <- as.numeric(substr(add$file,nchar(add$file),nchar(add$file)))
            add$file <- substr(add$file,1,(nchar(add$file)-1))
            df <- rbind(df,add)
          }
        }
      }
    }
    df <- df[which(!is.na(df$n)),]
  }else{   ### if age-independent dataset:
    for(i in c(grep('Total_Releases',names(k)), grep('Recovery_Data',names(k))))
    {
      print(names(k)[i])
      a0 <- k[[names(k)[i]]]
      a0[a0 == -1] <- NA
      if(grepl('Total_Releases',names(k)[i])){
        add <- data.frame(n=a0$n,y_rel=a0$year,y_rec=NA)
        add$file <- names(k)[i]
        add$age <- 1
        add$file <- substr(add$file,1,(nchar(add$file)-1))
        df <- rbind(df,add)
      }else{
        for(y in 1:nrow(a0)){
          for(x in 1:ncol(a0)){
            add <- data.frame(n=a0[y,x],y_rel=y_rel[y],y_rec=y_rec[x])
            add$file <- names(k)[i]
            add$age <- 1
            add$file <- substr(add$file,1,(nchar(add$file)-1))
            df <- rbind(df,add)
          }
        }
      }
    }
    if(age.dependent) warning('provided data set to remake is age independent! "age.dependent" was set to FALSE')
    age.dependent <- F
  }
  
  if(!any(!is.na(df$n[grep('Release_Tag_Recovery',df$file)])) & !only.harvested){
    warning('provided data set to remake does not contain tag released captures! "only.harvested" was set to TRUE')
    only.harvested <- T
  }
  # unique(df$file)
  
  
  #### create initial release file
  j <- grep('Total_Release',df$file)
  rel <- df[j,]
  rel$file <- rel$y_rec <- c()
  releases <- c()
  for(i in 1:nrow(rel)){
    rel.date <- as.Date(paste0(rel$y_rel[i],"-01-01"))-1+round(runif(rel$n[i],1,365))
    DeployID <- .create.serial.number(rel$n[i],8)
    add <- data.frame(DeployID=DeployID,rel.date=rel.date,age=rel$age[i])
    releases <- rbind(releases,add)
  }
  #   head(releases)
  
  
  recaptures <- c()
  rec <- df[-j,]
  rec <- rec[which(rec$n > 0),]
  rec$harvested <- grepl('Harvest',rec$file)
  rec$file <- c()
  for(i in 1:nrow(rec)){
    rel.date <- as.Date(paste0(rec$y_rel[i],"-01-01"))-1+round(runif(rec$n[i],1,365))
    rec.date <- as.Date(paste0(rec$y_rec[i],"-01-01"))-1+round(runif(rec$n[i],1,365))
    DeployID <- .create.serial.number(rec$n[i],8) # simplified, not linked releases-df
    add <- data.frame(DeployID=DeployID,rel.date=rel.date,rec.date=rec.date,age=rec$age[i],harvested=rec$harvested[i])
    recaptures <- rbind(recaptures,add)
  }
  
  releases$y_rel <- .as.year(releases$rel.date)
  recaptures$y_rel <- .as.year(recaptures$rel.date)
  recaptures$y_rec <- .as.year(recaptures$rec.date)
  
  #   inst.pkg(plyr)
  #   ddply(releases,.(y_rel),function(x)nrow(x))
  #   ddply(recaptures,.(y_rec),function(x)nrow(x))
  
  # head(releases)
  # head(recaptures)
  
  
  #   save(releases,recaptures,file="test.set_raw.rd")
  # save(total_rel_list,total_rec_list,only.harvested,y_rel,y_rec,ages,file="test.set.rd")
  #   load("test.set_raw.rd",verbose = T)
  
  #### create data input file:
  file.remove(paste0(new.setup.name,".dat"))
  make.dat(new.setup.name,releases,recaptures,age.dependent=age.dependent,only.harvested = only.harvested,...)
}



