clean.admb <- clean.IRATE <- function(setup.name){
  run_files <- Sys.glob(paste0(setup.name,'*'))
  i <- grep(paste0(setup.name,'.dat'),run_files)
  if(length(i) > 0) run_files <- run_files[-i]
  for(n in run_files) file.remove(n)
  file.remove('fmin.log')
  admb_files <- Sys.glob('admodel.*')
  for(n2 in admb_files) file.remove(n2)
  
}