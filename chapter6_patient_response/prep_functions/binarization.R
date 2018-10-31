#22.07.16 script for binarization of response value, 
# it uses approach from BHK paper (inflection point estimation when cor. coe is < 0.95) and some part of scripts

source("prep_functions/distancePointSegment.R")
source("prep_functions/distancePointLine.R")


binarization <- function(response_vect, grey_zone=F) {
  xx <- response_vect
  n <- which(is.na(xx))
  if (length(n)>0) {
    xx <- xx[-which(is.na(xx))]
  }
  if(grey_zone==F) {
    cor.min.linear=0.95
    oo <- order(xx, decreasing=TRUE)
    ## test linearity with Pearson correlation
    cc <- cor.test(-xx[oo], 1:length(oo), method="pearson")
    
    if(cc$estimate > cor.min.linear){
      ## approximately linear waterfall
      cutoff <- median(xx)
      
    } else {
      
      xx <- xx[-which(xx==max(xx))]
      xx <- xx[-which(xx==min(xx))]
      oo <- order(xx, decreasing=TRUE)
      ## line between the two extreme sensitivity values
      dd <- cbind("y"=xx[oo][c(1, length(oo))], "x"=c(1, length(oo)))
      rr <- lm(y ~ x, data=data.frame(dd))
      ## compute distance from sensitivity values and the line between the two extreme sensitivity values
      ddi <- apply(cbind(1:length(oo), xx[oo]), 1, function(x, slope, intercept) {
        return(distancePointLine(x=x[1], y=x[2], slope=slope, intercept=intercept))
      }, slope=rr$coefficients[2], intercept=rr$coefficients[1])
      ## non linear waterfall
      ## identify cutoff as the maximum distance
      cutoff <- (xx[oo])[which.max(abs(ddi))]
    }
    
    sens <- which(response_vect <= cutoff)
    res <- which(response_vect > cutoff)
    response_vect[sens] <- "sens"
    response_vect[res] <- "resist"
  } else {
    
    xx <- xx[-which(xx==max(xx))]
    xx <- xx[-which(xx==min(xx))]
    cutoffs <- quantile(xx, c(1/3,2/3))
    
    sens <- which(response_vect < cutoffs[1])
    res <- which(response_vect > cutoffs[2])
    interm <- setdiff(1:length(response_vect), c(sens,res))
    response_vect[sens] <- "sens"
    response_vect[res] <- "resist"
    response_vect[interm] <- NA
    
  }

return(response_vect)

}


####
