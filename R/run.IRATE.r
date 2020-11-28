run.IRATE <- function(setup.file,safe=F,re=F,verbose=T,
                      admb_errors=c("stop","warn","ignore"),
                      mcmc=F,mcmc.opts=mcmc.control(),
                      profile=F,extra.args=""){
  
  setup.file <- gsub('.dat',"",setup.file)
  file.copy(system.file("tplfiles/IRATEv2.tpl",package="IRATER"),paste0(setup.file,'.tpl')) # create tpl file
  
  R2admb::compile_admb(setup.file,safe=safe,re=re,verbose=verbose,admb_errors=admb_errors)
  
  R2admb::run_admb(setup.file,verbose=verbose,mcmc=mcmc,
           mcmc.opts=mcmc.opts,profile=profile,
           extra.args=extra.args,admb_errors=admb_errors)
}
