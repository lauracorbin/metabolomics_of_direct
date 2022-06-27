# functions for analysis script 04.0

# load libraries
library('lme4')
library('dplyr')
library("car")


# winsorise function
winsorize_x = function(x, cut=0.01) {
#winsorize_x = function(x, nsd=5) {
  cut_point_top <- quantile(x,1-cut,na.rm=T)
  cut_point_bottom <- quantile(x,cut,na.rm=T)
  #cut_point_top <- mean(x,na.rm=T) + (nsd*(sd(x,na.rm=T)))
  #cut_point_bottom <- mean(x,na.rm=T) - (nsd*(sd(x,na.rm=T)))
  i = which(x >= cut_point_top)
  x[i] = cut_point_top
  j = which(x <= cut_point_bottom)
  x[j] = cut_point_bottom
  return(x)
}

# convert to pa
pa.convert = function(dtst=dtst) {
  output <- apply(dtst,2,function(x) ifelse(is.na(x),0,1))
  metab_number <- (ncol(output)-10)/2
  for (j in 1:nrow(output)) {
    # set to NA if entire baseline sample missing
    if ( sum(output[j,1:metab_number] == 0) == metab_number ) {
      output[j,1:metab_number] <- NA
    }
  }
  return(as.data.frame(output))
}

