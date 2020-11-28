read.par <- function(fn)#,skip=0) 
{
  fn <- gsub('.par','',fn)
  fnext <- paste(fn,"par",sep=".")
  if (!file.exists(fnext)) stop("no PAR file found")
  
  tmp2 <- readLines(fnext)
  #   if (skip>0) tmp2 <- tmp2[-(1:skip)]
  tmp2[which(grepl('#\\*',tmp2))] <- " "
  parlines <- grep("^#",tmp2) # identify commented sections (long parameter names)
  vallines <- which(!grepl("^#",tmp2)) # identify commented sections (long parameter names)
  tmp2[parlines]
  npar <- length(parlines)  ## number of distinct parameters
  tmp2 <- gsub(':','',tmp2)
  
  parnames <- gsub("^# +","",gsub(":$","",tmp2[parlines])) # clean header
  
  p0 <- gsub("Number of parameters = ", "",parnames[1])
  p1 <- gsub(" Objective function value = ", "#",p0)
  p2 <- gsub("  Maximum gradient component = ", "#",p1)
  add.vals <- as.numeric(strsplit(p2,'#')[[1]])
  
  vals.list <- list()
  vals.list$nparam <- add.vals[1]
  vals.list$ofv <- add.vals[2]
  vals.list$max_grad_comp <- add.vals[3]
  
  parlines2 <- c(parlines,length(tmp2))
  for(i in 2:length(parlines)){
    j <- (parlines2[i]+1):(parlines2[i+1]-1)
    j <- j[which(j > (parlines2[i]))]
    vals0 <- c()
    for(jj in j) {
      add0 <- strsplit(tmp2[jj]," ")
      add <- as.numeric(add0[[1]][which(nchar(add0[[1]]) > 0)])
      vals0 <- rbind(vals0,add)
    }
    rownames(vals0) <- c()
    vals.list[[parnames[i]]] <- vals0
  }

  return(vals.list)
}