#' Calculate group-specific covariance and correlation matrices
#'
#' @param X the data matrix
#' @param y the label
#'
#' @return COR.list is a list, each element contains a group-specific correlation matrix
#' @return COV.list is a list, each element contains a group-specific covariance matrix
#' @return COV.scale.list is a list, each element contains a standardized group-specific covariance matrix
#' @export
individual_covcor <- function(X ,y ){
  K = length(unique(y))
  marker.name.vec <- colnames(X )
  p = dim(X)[2]
  if(is.null(marker.name.vec)){
    marker.name.vec <- c(1:p) #paste("V", c(1:p), sep="")
    colnames(X) <- marker.name.vec
  }


  COV.list <- list()
  for(k in 1:K){print(k)
    idk<- which(y==k);
    COV.tmp  <- stats::cov(X[idk, ])
    colnames(COV.tmp) <- marker.name.vec
    rownames(COV.tmp) <- marker.name.vec
    COV.list[[k]]   <- COV.tmp
  }

  COV.scale.list <- list()
  for(k in 1:K){print(k)
    idk<- which(y==k);
    COV.scale.tmp  <- stats::cov(scale(X[idk, ]))
    colnames(COV.scale.tmp) <- marker.name.vec
    rownames(COV.scale.tmp) <- marker.name.vec
    COV.scale.list[[k]]   <- COV.scale.tmp
  }

  COR.list <- list()
  for(k in 1:K){ print(k)
      idk<- which(y==k);
      COR.tmp  <- stats::cor(X[idk, ])
      colnames(COR.tmp) <- marker.name.vec
      rownames(COR.tmp) <- marker.name.vec
      COR.list[[k]]   <- COR.tmp
    }

  for(k in 1:K){ print(k)
      COR.tmp  <- COR.list[[k]]
      COR.tmp[which(is.na(COR.tmp), arr.ind=TRUE)] <- 0
      COR.list[[k]]   <- COR.tmp
  }



  return(list(COR.list = COR.list, COV.list = COV.list,
              COV.scale.list = COV.scale.list))
}
