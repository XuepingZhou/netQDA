ConnectComponent <-  function( Screenset  ,
                               COR ,  d , nb ){
  p = dim(COR)[1]
  marker.name.vec <- colnames(COR)
  vis=rep(0, p)
  id <- Screenset[1]
  scr <- Screenset

  local.network_list <- list()

  while(!is.na(id)){
    tmp <- con_node_dn(id, COR, visited=vis, depth.max=d, neighbor.max=nb)
    con_set <- tmp$connect
    names(con_set) <- marker.name.vec[con_set]
    vis <- tmp$visited
    local.network_list <- append(local.network_list, list(con_set))
    id <- scr[is.na(match(scr, con_set))][1]
    scr <- scr[is.na(match(scr, con_set))]
  }
  return(local.network_list)
}
