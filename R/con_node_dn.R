con_node_dn <-
  function(vid, Sig, visited=rep(0, dim(Sig)[1]), depth=0, depth.max=2, neighbor.max=10, connect=NULL, path=NULL, depth.path=NULL){

    if(dim(Sig)[1] != dim(Sig)[2]){ stop("Sig is not a square matrix") }
    if(vid>dim(Sig)[1] | vid<=0){ stop("Wrong vertex index")}

    neighbor <- 0

    if(neighbor<=neighbor.max & depth<=depth.max){
      connect <- c(connect, vid)
      path <- c(path, vid)
      visited[vid] <- 1
    }


    flag <- 0

    tmp <- which(visited==0)
    ranks <- order(-abs(Sig[vid,tmp]))
    rest.ranks <- tmp[ranks]

    depth0 <- depth
    neighbor0 <- neighbor


    for(i in rest.ranks){
      depth1 <- depth0

      if(Sig[vid,i]!=0 & visited[i]==0 & neighbor<neighbor.max & depth1<depth.max){



        neighbor <- neighbor +1


        depth1 <- depth1+1

        tmp <- con_node_dn(i, Sig, visited=visited, depth=depth1, depth.max=depth.max, neighbor.max=neighbor.max, connect=connect, path=path, depth.path=depth.path)

        visited <- tmp$visited
        depth <- tmp$depth
        connect <- tmp$connect
        path <- tmp$path

        depth.path <- tmp$depth.path
        depth.path <- c(depth1, depth.path)

        flag <-1

      }  # end of  if(Sig[vid,i]!=0 & visited[i]==0 & neighbor<=neighbor.max)
    } # end of for



    if(flag==0){stop}

    return(list(connect=connect, visited=visited, depth=depth, path=path, depth.path=depth.path))
  }
