Omega_sub_cal <- function(
    group.list , COV
){

  Omega.s <- matrix(0, dim(group.list)[1], dim(group.list)[1])
  COV.s <- matrix(0, dim(group.list)[1], dim(group.list)[1])
  colnames(Omega.s) <- group.list$component
  rownames(Omega.s) <- group.list$component
  colnames(COV.s) <- group.list$component
  rownames(COV.s) <- group.list$component


  det = 1
  det.all <- c()
  for(i.group in unique(group.list$group)){
    con_set <- group.list$component[group.list$group == i.group]
    COV.s[as.character(con_set), as.character(con_set)] <- COV[con_set, con_set]
    Omega.s[as.character(con_set), as.character(con_set)] <- my.inv(COV[con_set, con_set])

    tmp.det<- det(as.matrix(COV[con_set, con_set]))
    det.all <- c(det.all, tmp.det)
    if(tmp.det > 0 ){

      det = tmp.det * det
    }

  }

  return(list(Omega.s = Omega.s, det.cov = det,  det.all =  det.all))
}

