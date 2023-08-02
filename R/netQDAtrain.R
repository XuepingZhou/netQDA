#' netQDA Train Step
#'
#' @param X is a n by p matrix
#' @param y is a vector of length n
#' @param COV.list is a list with each element corresponding to the covariance matrix of each group.
#' @param COR.list is a list with each element corresponding to the correlation matrix of each group.
#' @param COV.scale.list is a list with each element corresponding to the covariance matrix of standadized data in each group.
#' @param pair the index for group. For example, pair=c(1,2), denotes the comparison between group 1 and 2
#' @param tau is the hyper-parameter for the number of strong differential expressed genes.
#' @param alpha is the hyper-parameter for the threshold for correlation matrix. We suggest using threshold no smaller than 0.7
#' @param nu is the hyper-parameter for the number of differential expressed (DE) genes.
#' @param d is the hyper-parameter depth/layer when searching for connected component based on the strong differential expressed (DE) gene.
#' @param nb is the hyper-parameter of number of features for each layer when searching for connected component based on the strong differential expressed gene.
#' @param delta is the hyper-parameter for the number of differential connected (DC)  genes.
#' @param p1 is the probability of a data point belonging to group 1.
#' @param p2 is the probability of a data point belonging to group 2. p1 + p2 = 1.
#'
#' @return try.meanDiff.firstTau the mean differences of data in group 1 and group 2
#' @return try.local.DCnetwork_list connected component of DC features based on DE features
#' @return try.local.DEnetwork_list connected component of DE features based on DE features
#' @return screenset.DC selected DC features
#' @return screenset.DE selected DE feautres
#' @return screenset.DE.iff if and only if condition values for DE features
#' @return screenset selected features (both DC and DE)
#' @return Omega.cis.mat.list Precision matrix of each connected component of selected features.
#' @return det.cov.mat determinants of covariance matrix of each group.
#' @return log.det.cov.mat log of determinants of covariance matrix of each group.
#' @return p.vec the probability of belong to group 1 and group 2.
#' @return pair the group index of the two groups.
#' @importFrom stats
#' @export
netQDAtrain <- function(X, y ,
                        COV.list=NULL,  COR.list=NULL,    COV.scale.list=NULL,
                        pair = c(1,2),
                        tau=20, alpha=0.7,  nu= 50,
                        d=2, nb=10,
                        delta = 100  ,
                        p1 = NULL, p2 = NULL

){
 K =2 # binary classification
  ## ID ---------------------------------------
  p <- dim(X)[2]
  n <- dim(X)[1]

  if(length(y)!=n){ stop("X and Y contain different number of samples!!") }

  id1<- which(y==pair[1])
  if(length(id1)==0){
    stop(paste("There is no y entries labeled as ", pair[1], "!", sep=""))	}

  id2<- which(y==pair[2])
  if(length(id2)==0){
    stop(paste("There is no y entries labeled as ", pair[2], "!", sep=""))}

  marker.name.vec <- colnames(X)
  if(is.null(marker.name.vec)){
    marker.name.vec <- c(1:p) #paste("V", c(1:p), sep="")
    colnames(X) <- marker.name.vec
  }

  y12 <- c(rep(0, length(id1)), rep(1, length(id2)))
  X12 <- X[c(id1, id2),]

  if( sum(duplicated(colnames(X))) > 0) {
    stop(paste("X contains duplicated colnames!"))}

  if(is.null(COR.list) |  is.null(COV.list) | is.null(COV.scale.list)){
    tmp <- individual_covcor(X=X, y=y)
    COR.list <-tmp[["COR.list"]]
    COV.list <-tmp[["COV.list"]]
    COV.scale.list <-tmp[["COV.scale.list"]]
  }else{
    COR1 <- COR.list[[pair[1]]][marker.name.vec , marker.name.vec ]
    COR2 <- COR.list[[pair[2]]][marker.name.vec , marker.name.vec ]
    COV1 <- COV.list[[pair[1]]][marker.name.vec , marker.name.vec ]
    COV2 <- COV.list[[pair[2]]][marker.name.vec , marker.name.vec ]
    COV1.scale <- COV.scale.list[[pair[1]]][marker.name.vec , marker.name.vec ]
    COV2.scale <- COV.scale.list[[pair[2]]][marker.name.vec , marker.name.vec ]
  }



  ## DE  ---------------------------------
  mu1 = apply(X12[which(y12==0),], 2, mean)
  mu2 = apply(X12[which(y12==1),], 2, mean)
  meanDiff <- mu1 - mu2
  meanDiffRank <- order(-abs(meanDiff))
  meanDiff_firstTau <- meanDiffRank[1:tau];
  names(meanDiff_firstTau) <- marker.name.vec[meanDiff_firstTau]
  meanAvg <- (mu1 + mu2) /2
  meanDiff[meanDiff_firstTau]


  ## DE-connected features ---------------------------------
  COR1[which(abs(COR1)<alpha, arr.ind=TRUE)] <- 0
  COR2[which(abs(COR2)<alpha, arr.ind=TRUE)] <- 0
  local.DEnetwork_list1 <- ConnectComponent(Screenset =meanDiff_firstTau  ,
                                            COR = COR1,  d = d, nb = nb)
  local.DEnetwork_list2 <- ConnectComponent(Screenset =meanDiff_firstTau  ,
                                            COR = COR2,  d = d, nb = nb)

  local.DEnetwork_list <- c( local.DEnetwork_list1 ,  local.DEnetwork_list2)
  for(a in 1:(length(local.DEnetwork_list) )){
    local.DEnetwork_list[[a]] <- names(local.DEnetwork_list[[a]])
  }
  combine.group <- UniqueComponent( local.DEnetwork_list, screenset = NULL)

  Omega1.s <-  Omega_sub_cal(group.list =  combine.group , COV = COV1)[[1]]
  Omega2.s <-  Omega_sub_cal(group.list =  combine.group , COV = COV2)[[1]]

  Omega1.scale.s <-  Omega_sub_cal(group.list =  combine.group , COV = COV1.scale)[[1]]
  Omega2.scale.s <-  Omega_sub_cal(group.list =  combine.group , COV = COV2.scale)[[1]]


  # JI feature ------------------------------------
  iffcond.s <-   as.vector(Omega2.s%*% meanDiff[colnames(Omega2.s)])
  names(iffcond.s) <-colnames(Omega2.s)
  # iffcond.s <- iffcond.s[abs(iffcond.s) > 1e-5] ###### //
  iffcond.s_order <- iffcond.s[order(abs(iffcond.s), decreasing = T)]

  if(length(iffcond.s_order)>nu){
    screenID.DE <- names(iffcond.s_order)[1:nu]
    iffcond.s_order_topnu = iffcond.s_order[1:nu]
  }else{
    screenID.DE <- names(iffcond.s_order)
    iffcond.s_order_topnu = iffcond.s_order
  }

  ## DC features ---------------------------------------
  omega.diff.delta = NULL
  local.DCnetwork_list <-  NULL
  screenID.DC = NULL

  if(delta > 0 ){
    Omega.diff.half= abs(Omega1.scale.s - Omega2.scale.s)
    Omega.diff.half[lower.tri(Omega.diff.half , diag = T)] <- 0
    delta = min(c(delta, dim(Omega.diff.half)[2]))
    omega.diff.delta <- sort(c(abs( Omega.diff.half )),
                             decreasing = T)[1:delta]

    if(length(omega.diff.delta) >0 ){
      local.DCnetwork_list <- list()
      if(sum( omega.diff.delta > 1e-2 )  > 0 ){
        for( dd in 1:sum( omega.diff.delta > 1e-2 )  ){ #####/////////
          i.pos.tmp <-which( Omega.diff.half ==  omega.diff.delta[dd], arr.ind = TRUE)
          local.DCnetwork_list[[dd]] <- colnames( Omega.diff.half)[as.vector(i.pos.tmp)]
        }
      }
      screenID.DC<- sort(unique( unlist( local.DCnetwork_list)))
    }
  }

  ## Final screenset ----------------------------------------------
  screenset = unique(sort(c(screenID.DC, screenID.DE)))  # print(screenset)
  local.network_list <-  append(local.DCnetwork_list,
                                local.DEnetwork_list)
  combine.group <- UniqueComponent( local.network_list,
                                    screenset =   screenset)
  Omega1.s <-  Omega_sub_cal(group.list =  combine.group , COV = COV1)[[1]]
  Omega2.s <-  Omega_sub_cal(group.list =  combine.group , COV = COV2)[[1]]

  ## summary
  mean.cis.mat <- NULL
  for(k in 1:K){
    id <- which(y==k)
    if(length(screenset) > 1){
      X.k <- X[id,  screenset,drop = F]
      mean.cis.vec <- apply(X.k, 2, mean)
    }else{
      X.k <- X[id,  screenset]
      mean.cis.vec <- mean(X.k )
    }
    mean.cis.mat <- cbind(mean.cis.mat, mean.cis.vec)
  }

  Omega.cis.mat.list <- list()
  log.det.cov.mat <- c()
  det.cov.mat <- c()
  for(k in 1:K){
    COVk = COV.list[[k]]
    tmp.omega.det <-  Omega_sub_cal(group.list =  combine.group , COV = COVk)
    Omega.cis.mat <-  tmp.omega.det[[1]]
    Omega.cis.mat.list[[k]] <- Omega.cis.mat
    tmp <-  tmp.omega.det[[2]]
    if(tmp >0 ){
      det.cov.mat <- c(  det.cov.mat,  tmp )
      log.det.cov.mat <- c( log.det.cov.mat ,  log( det.cov.mat[k] ))
    }else{
      det.cov.mat <- c(  det.cov.mat, 0)
      log.det.cov.mat <- c( log.det.cov.mat , -9999)
    }
  }

  if(is.null(p1) & is.null(p2)){
    p1 = length(id1)/(length(id1)+length(id2))
    p2=  length(id2)/(length(id1)+length(id2))
  }

  # Return
  result <- list(
    try.meanDiff.firstTau = meanDiff[meanDiff_firstTau],
    try.local.DCnetwork_list= local.DCnetwork_list,
    try.local.DEnetwork_list= local.DEnetwork_list,

    screenset.DC = screenID.DC,
    screenset.DE = screenID.DE,
    screenset.DE.iff =   iffcond.s_order_topnu ,
    screenset=screenset ,
    mean.cis.mat = mean.cis.mat ,
    Omega.cis.mat.list = Omega.cis.mat.list ,
    det.cov.mat =  det.cov.mat ,
    # det.cov.mat.all =  tmp.omega.det[[3]],
    log.det.cov.mat =  log.det.cov.mat ,
    p.vec =  c(p1,p2),    pair = pair
  )
  result
}
