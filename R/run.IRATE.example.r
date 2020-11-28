run.IRATE.example <- function(setup.name, run.folder,...){
  
  owd <- getwd()
  
  ### select setup:
  setups <- IRATE.examples()
  valid.setup <- F
  if(!missing(setup.name)) {
    valid.setup <- setup.name %in% setups # check if user defined setup is valid
    if(!valid.setup) warning("user defined setup '",setup.name,"' is not valid!")
  }
  if(missing(setup.name) | !valid.setup){ # ask to select valid setup
    i <- readline(paste("\nPlease choose one of the following example setups to run by selecting its number (see: for documentation), \n",
                        paste(paste0("[",1:length(setups),"]"),setups,collapse=", "),"\n"))
    i <- as.numeric(i)
    fn <- setups[i]
  }else{
    fn <- setup.name
  }
  if(!(fn %in% setups)) stop("User defined setup is not valid! Please rerun function and check ?run.IRATE.example or IRATE.examples() for valid IRATE example setups")
  
  ### create run folder
  if(missing(run.folder)) run.folder <- paste0(getwd(),'/',fn)
  cat('creating run.folder:', run.folder)
  dir.create(run.folder,recursive = T)
  dat.file <- system.file(paste0("IRATE.examples/",fn,".dat"),package="IRATER")
  file.copy(dat.file, run.folder)
  
  setwd(run.folder)
  try(run.IRATE(fn))
  setwd(owd)
}

IRATE.examples <- function(){
  owd <- getwd()
  examples.folder <- paste0(system.file("IRATE.examples",package="IRATER"),"/")
  setwd(examples.folder)
  setups <- gsub('*.dat','',Sys.glob('*.dat'))
  setwd(owd)
  return(setups)
}