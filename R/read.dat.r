read.dat <- function(setup.file.name, skip=0) {
  setup.file.name <- gsub('.dat','',setup.file.name)
  fnext <- paste(setup.file.name,"dat",sep=".")
  if (!file.exists(fnext)) stop("no Input-file found")
  
  tmp2 <- readLines(fnext)
  if (skip>0) tmp2 <- tmp2[-(1:skip)]
  tmp2[which(grepl('#\\*',tmp2))] <- " "
  parlines <- grep("^#",tmp2) # identify commented sections (long parameter names)
  tmp2[parlines]
  npar <- length(parlines)  ## number of distinct parameters
  
  pp <- c(parlines,length(tmp2)+1)
  istart <- which(diff(pp) > 1)[1]
  
  parlist <- list()
  parlist$data_info <- as.matrix(tmp2[1:(istart-1)],ncol=1)
  
  parnames0 <- gsub("^# +","",gsub(":$","",tmp2[parlines])) # clean header
  parnames <- gsub(' ','_',parnames0)
  parnames <- gsub('-','_',parnames)
  parnames <- gsub('_\\&_',"_", parnames)
  parnames <- gsub('_\\/_',"_", parnames)
  parnames <- gsub('___',"_", parnames)
  parnames <- gsub('__',"_", parnames)
  parnames <- gsub('Number_of_ages','Number_of_Ages',parnames)
  
  for(i in 1:length(parnames)) {
    pn <- parnames[i]
    if(substr(pn,nchar(pn),nchar(pn)) == "_") {
      parnames[i] <- substr(pn,1,(nchar(pn)-1))
      i <- i - 1
    }
  }
  
  yrel <- yrec <- c()
  
  for (ii in istart:npar) {
    nrows <- diff(pp)[ii]-1 # of variable entries
    if(nrows > 0){
      i <- parlines[ii]
      
      do.true.false <- ((grepl('Combine', parnames[[ii]]) | grepl('Estimate', parnames[[ii]]) | grepl('\\?', parnames[[ii]])) | grepl('Use', parnames[[ii]]))
      
      if(grepl('Estimate', parnames[[ii]]) & grepl('\\?', parnames[[ii]])){
        parnames[[ii]] <- strsplit(parnames[[ii]],"\\?")[[1]][1]
      }
      
      parlist[[parnames[ii]]] <- as.matrix(tmp2[(i+1):(i+nrows)],ncol=1)
      
      s <- 1
      add <- c()
      for(n in 1:nrows){
        add0 <- tmp2[(i+s):(i+n)]
        add0 <- gsub('  ',' ',add0)
        add0 <- as.numeric(strsplit(add0,' ')[[1]])
        if(do.true.false)  add0 <- .onekey2truefalse(add0)
        if(any(!is.na(add0))){
          if(nrows > 1) {
            add <- rbind(add,add0)
            rownames(add) <- c()
          }else{
            add <- add0
          }
          s <- s+1
        }
      }
      if(!is.null(dim(add))) {
        if(all(is.na(add[,1]))) add <- add[,-1]
      }
      if(grepl("Total_Releases",parnames[ii])){
        yrel <- data.frame(year=parlist$Release_Start_End_Years[1]:parlist$Release_Start_End_Years[2])
        yrec <- data.frame(year=parlist$Recapture_Start_End_Years[1]:parlist$Recapture_Start_End_Years[2])
        add <- data.frame(year=yrel,n=add)
      }
      if(parnames[ii] %in% c("F_Starting_Values")){
        y <- data.frame(year=parlist$F_Period_Beginning_Year)
        add <- data.frame(year=y$year,F=add)
      }
      if(parnames[ii] %in% c("Hooking_Mortality_Rates")){
        add <- data.frame(year=yrec,F=add)
      }
      if(parnames[ii] %in% c("Nonmixing_Starting_Values")){
        y <- data.frame(year=parlist$NM_Periods_Beginning_Years)
        add <- data.frame(year=y$year,F=add)
      }
      if(!is.null(dim(add))) {
#         print(dim(add))
#         print(length(yrec$year))
        if(!any(!(dim(add) == c(length(yrel$year),length(yrec$year))))){
          rownames(add) <- paste0("rel_",yrel$year)
          colnames(add) <- paste0("rec_",yrec$year)
        }
      }
      
      parlist[[parnames[ii]]] <- add
    }
  }
  
  return(parlist)
}

.onekey2truefalse <- function(x){
  if(x == -1) x <- F
  if(x == 1) x <- T
  return(x)
}

