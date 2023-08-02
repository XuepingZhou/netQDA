#' netQDA Classification Step
#'
#' @param TrainModel the output from netQDAtrain function.
#' @param X.new new data for classification task.
#' @return  Fisher.matrix Fisher's matrix
#' @return  predClass predicted group membership
#' @return  p1 predicted probability of group 1
#' @return  p2 predicted probability of group 2
#' @export
netQDApredict <- function(TrainModel , X.new  ){

  mean.cis.mat = TrainModel[["mean.cis.mat"]]
  Omega.cis.mat.list = TrainModel[["Omega.cis.mat.list"]]
  screenset = TrainModel[["screenset"]]
  det.cov.mat = TrainModel[["det.cov.mat"]]
  log.det.cov.mat = TrainModel[["log.det.cov.mat"]]
  pair = TrainModel[["pair"]]
  p.vec <-  TrainModel[["p.vec"]]
  logp1 = log( p.vec[1]); logp2 = log( p.vec[2])

  p = length(screenset)
  n = dim(X.new)[1]

  mu1 = as.matrix(mean.cis.mat[, 1] , nrow=1, ncol=length(screenset))
  mu2 = as.matrix(mean.cis.mat[, 2] , nrow=1, ncol=length(screenset))

  mu1.matrix = matrix(rep(mu1,n), nrow=n, ncol = length(screenset ), byrow = T )
  mu2.matrix = matrix(rep(mu2,n), nrow=n, ncol = length(screenset ), byrow = T )

  offset1 <- matrix(-0.5*  log.det.cov.mat[1]  + logp1 , ncol = 1, nrow = n)
  offset2 <- matrix(-0.5*  log.det.cov.mat[2]  + logp2,  ncol = 1, nrow = n)

  DEV1 <-  diag(0.5* as.matrix(X.new[, screenset, drop= F] - mu1.matrix) %*%
                  Omega.cis.mat.list[[1]]  %*%
                  t(as.matrix(X.new[, screenset, drop= F] - mu1.matrix)) )
  DEV2 <-  diag(0.5* as.matrix(X.new[, screenset, drop= F] - mu2.matrix) %*%
                  Omega.cis.mat.list[[2]]  %*%
                  t(as.matrix(X.new[, screenset, drop= F] - mu2.matrix)))

  ## ADD end

  Fisher.vec1 <-   offset1  -   DEV1 -  length(screenset) /2 * log(2*pi)
  Fisher.vec2 <-   offset2  -   DEV2 -  length(screenset) /2 * log(2*pi)

  Fisher.matrix = cbind(Fisher.vec1, Fisher.vec2)
  predClass <- apply(Fisher.matrix, 1, function(x)  if(x[1] > x[2]){pair[1]}else{pair[2]} )
  # data <- softmax(data = Fisher.matrix)

  ##### posterior
  p1 <- fx(X = X.new[, screenset, drop=F],  DEV= DEV1, detcovk=  det.cov.mat[1])
  p2 <- fx(X = X.new[, screenset, drop=F],  DEV= DEV2, detcovk=  det.cov.mat[2] )

  p1.true = (p1 * p.vec[1] ) / (p1 * p.vec[1] +  p2 * p.vec[2] )
  p2.true = (p2 * p.vec[2] ) / (p1 * p.vec[1] +  p2 * p.vec[2] )

  for(ii in which(is.na(p2.true))){
    if(Fisher.vec2[ii] > Fisher.vec1[ii]){
      p2.true[ii] = 1
      data$x2[ii] = 1
    }else{
      p2.true[ii] = 0
      data$x2[ii] = 0
    }

  }

  result <- list( Fisher.matrix=Fisher.matrix,
                  predClass=predClass,
                  # probClass = data$x2,
                  p1=p1.true, p2=p2.true )

}
