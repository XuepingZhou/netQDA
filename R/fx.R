fx <- function(X, DEV, detcovk){
  p = dim(X)[2] ; n = dim(X)
  tmp = 1/ ((2*pi)^(p/2) * sqrt(detcovk)) * exp(-1* DEV)
  return( tmp)
}
