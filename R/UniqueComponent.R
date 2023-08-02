UniqueComponent <- function(componentlist ,
                            screenset = NULL ){

  for(a in 1:(length(componentlist) )){
    avec <- componentlist[[a]]
    if(!is.null(screenset )){
      avec <-   avec [  avec  %in%  screenset ==T] ###
    }
    componentlist[[a]] <- avec
  }



  componentlist <-  componentlist[lengths(componentlist) > 0L]
  if(length(unlist(componentlist)) ==
     length(unique( unlist(componentlist)))){
    component  <- c() ; group <- c()
    for(i in 1:length( componentlist)){
      component_current <- componentlist[[i]]

      group_current <- rep(i, length(component_current))
      component <- c(component , component_current) ;
      group <- c(group, group_current)
    }
    component.group.df <- as.data.frame(cbind(component = component,
                                              group= group))
    return(component.group.df)

  }else{
    B = length(componentlist)
    for(a in 1:(B -1)){
      avec <- componentlist[[a]]

      if ( length(avec) >0 ){
        for( b in  (a+1):B ) {
          bvec <- componentlist[[b]]
          if(length(intersect(avec, bvec)) !=0 ){
            avec = unique(c(avec, bvec));
            componentlist[[b]] <- character(0)
          }
        }
      }
      componentlist[[a]] <- sort(avec)
    }
    componentlist <-  componentlist[lengths(componentlist) > 0L]
    UniqueComponent(componentlist)
  }
}

