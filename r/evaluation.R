### distance -------------------------------------------------------------------------------
Dis <- function(Z1,Z2){
  d = base::norm(Z1%*%t(Z1)-Z2%*%t(Z2),type='F')^2 / nrow(Z1)^2
  return(d)
}


mDis <- function(Z1,Z2){
  Ts = dim(Z1)[3]
  dis = rep(NA,Ts)
  for (t in 1:Ts) {
    if(length(dim(Z2))==3){
      dis[t] = base::norm(Z1[,,t]%*%t(Z1[,,t])-Z2[,,t]%*%t(Z2[,,t]),type='F')^2
    }else if(length(dim(Z2))==2){
      dis[t] = base::norm(Z1[,,t]%*%t(Z1[,,t])-Z2%*%t(Z2),type='F')^2
    }
  }
  return(mean(dis)/dim(Z1)[1]^2)
}

