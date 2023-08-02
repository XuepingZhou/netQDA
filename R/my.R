my.inv <-
  function(X, eps=1e-12){
    eig.X <- eigen(X, symmetric=TRUE)
    P <- eig.X[[2]]
    lambda <- eig.X[[1]]
    ind <- lambda > eps
    lambda[ind] <- 1/lambda[ind]
    lambda[!ind] <- 0
    ans <- P%*%diag(lambda, nrow=length(lambda))%*%t(P)
    return(ans)
  }
