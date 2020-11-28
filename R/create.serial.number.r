
.create.serial.number <- function(n,digits){
  out <- c()
  for(x in 1:n){
    add <- c()
    for(i in 1:digits){
      j <- round(runif(1,1,3))
      add0 <- switch(j,
                     as.character(round(runif(1,1,9))),
                     LETTERS[round(runif(1,1,26))],
                     letters[round(runif(1,1,26))]
      )
      add <- paste0(add,add0)
    }
    out <- c(out,add)
  }
  return(out)
}