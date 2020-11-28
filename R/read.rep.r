read.short.rep <- function(fn,skip=0,show.sm=F,short=T) {read.rep(fn,skip=skip,show.sm=show.sm,short=short)}
  
read.rep <- function(fn,skip=0,show.sm=F,short=F) {
  fn <- gsub('.rep','',fn)
  fnext <- paste(fn,"rep",sep=".")
  if (!file.exists(fnext)) stop("no REP file found")
  
  dat <- read.dat(fn)
  yrec <- data.frame(year=dat$Recapture_Start_End_Years[1]:dat$Recapture_Start_End_Years[2])
  yrel <- data.frame(year=dat$Release_Start_End_Years[1]:dat$Release_Start_End_Years[2])
  
  tmp2 <- readLines(fnext)
  if (skip>0) tmp2 <- tmp2[-(1:skip)]
  tmp2[which(grepl('#\\*',tmp2))] <- " "
  parlines <- grep("^#",tmp2) # identify commented sections (long parameter names)
  tmp2[parlines]
  npar <- length(parlines)  ## number of distinct parameters
  
  pp <- c(parlines,length(tmp2)+1)
  istart <- which(diff(pp) > 1)[1]
  
  parlist <- list()
  parlist$rep_info <- as.matrix(tmp2[1:(istart-1)],ncol=1)
  
  parnames0 <- gsub("^# +","",gsub(":$","",tmp2[parlines])) # clean header
  parnames <- gsub(' ','_',parnames0)
  parnames <- gsub('-','_',parnames)
  for(p in 1:length(parnames)) {
    pn <- parnames[p]
    if(substr(pn,nchar(pn),nchar(pn)) == "_") parnames[p] <- substr(pn,1,(nchar(pn)-1))
  }
  for (ii in istart:npar) {
    nrows <- diff(pp)[ii]-1 # of variable entries
    i <- parlines[ii]
    
    parlist[[parnames[ii]]] <- as.matrix(tmp2[(i+1):(i+nrows)],ncol=1)
    
    s <- 1
    add <- c()
    for(n in 1:nrows){
      add0 <- tmp2[(i+s):(i+n)]
      add0 <- as.numeric(strsplit(add0,' ')[[1]])
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
    if(!is.null(dim(add))) if(all(is.na(add[,1]))) add <- add[,-1]
    if(!is.null(dim(add))) {  
      if(!any(!(dim(add) == c(length(yrel$year),length(yrec$year))))){
        rownames(add) <- paste0("rel_",yrel$year)
        colnames(add) <- paste0("rec_",yrec$year)
      }
    }
    if(parnames[ii] == "Total_Released"){
      add <- data.frame(year=yrel,n=add)
    }
#     if(parnames[ii] == "Total_Recovered_Tags"){
#       add <- data.frame(year=yrec,n=add)
#     }
    
    parlist[[parnames[ii]]] <- add
  }
  
  
  output_short <- list()
  stats <- list(c('Log_L',"AIC",'AICc','Effective_Sample_Size'),
                c('Unpooled_Chi_square','Upooled_df','Unpooled_c_hat'),
                c('Pooled_Chi_square','Pooled_df','Pooled_c_hat'))
  for(i in 1:length(stats)){
    sm <- list()
    for(field in stats[[i]]) sm[[field]] <- output_short[[field]] <- parlist[[field]]
    sm <- as.data.frame(sm)
    rownames(sm) <- ""
    if(show.sm) print(sm)
  }
  
  if(short) parlist <- output_short
  return(parlist)
}