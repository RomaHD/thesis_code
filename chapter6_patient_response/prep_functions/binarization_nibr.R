# 1.08.16
# very simple binarization func. for nibr
binarization_nibr <- function(response_vect, grey_zone=F) {

  cutoff=0
  sens <- which(response_vect < cutoff)
  res <- which(response_vect > cutoff)
  interm <- which(response_vect == cutoff)
  response_vect[sens] <- "sens"
  response_vect[res] <- "resist"
  response_vect[interm] <- "resist"
  
  if (grey_zone==T)
  {response_vect[interm] <- "NA"}
  return(response_vect)
  
}