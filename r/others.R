pca_est <- function(Y, q){
  p = ncol(Y)
  Ts = nrow(Y)
  Ft = sqrt(Ts)*as.matrix(eigen(Y%*%t(Y))$vectors)[,1:q]
  B = Ts^(-1) * t(Y)%*%Ft
  BF = Ft%*%t(B)
  return(list(pre_F=Ft, pre_B=B, pre_BF=BF)) 
}



local_pca <- function(Y, q){
  p = ncol(Y)
  Ts = nrow(Y)
  Ys = array(data=NA, dim=c(Ts,p,Ts))
  B_tilde = array(data=NA, dim=c(p,q,Ts))
  pre_F = matrix(data=NA, nrow=Ts, ncol=q)
  BF = matrix(data=NA, nrow=Ts, ncol=p)
  for (s in 1:Ts){
    for (t in 1:Ts){
      Ys[t,,s] = Y[t,]*k(Ts,p,t,s)^(1/2)
    }
    Fs_tilde = sqrt(Ts)*as.matrix(eigen(Ys[,,s]%*%t(Ys[,,s]))$vectors)[,1:q]
    B_tilde[,,s] = Ts^(-1) * t(Ys[,,s])%*%Fs_tilde
    if(s>1){
      for(j in 1:q){
        if(sign(B_tilde[1,j,1])!=sign(B_tilde[1,j,s]))
          B_tilde[,j,s] = B_tilde[,j,s] * (-1)
      }
    }
    pre_F[s,] = solve(t(B_tilde[,,s])%*%B_tilde[,,s]) %*% t(B_tilde[,,s])%*%Y[s,]
    BF[s,] = B_tilde[,,s]%*%pre_F[s,]
  }
  return(list(B_tilde=B_tilde, pre_F=pre_F, BF=BF))
}


TV_WPLL <- function(Y, Z=NULL, q){
  # Y: T by P matrix
  Y <- as.matrix(Y)
  Ts = nrow(Y)
  p = ncol(Y)
  #h = 1.059 * Ts^(-1/5) 
  h = (2.35/sqrt(12)) * Ts^(-1/5) * p^(-1/10)
  if(is.null(Z)){
    Z=matrix((1:Ts)/Ts,nrow=Ts)
  }
  Z = as.matrix(Z)
  Ys = array(data=NA, dim=c(Ts, p, Ts))
  B = array(data=NA,dim=c(p, q, Ts))
  Fac = array(data=NA,dim=c(Ts, q))
  BF = array(data=NA,dim=c(Ts, p))
  for (t in 1:Ts) {
    k = array(NA,Ts)
    for (s in 1:Ts) {
      k[s] = h^{-1} * Kernel(sqrt(sum((Z[s,]-Z[t,])^2))/h)
    }
    k = k/mean(k)
    for (s in 1:Ts){
      Ys[s,,t] = Y[s,]*k[s]^(1/2)
    }
    Fs_tilde = sqrt(Ts)*as.matrix(eigen(Ys[,,t]%*%t(Ys[, , t]))$vectors)[,1:q]
    B[, , t] = Ts^(-1) * t(Ys[,,t])%*%Fs_tilde
    Fac[t,] = solve(t(B[, , t])%*%B[, , t]) %*% t(B[, , t])%*%Y[t,]
    BF[t,] = B[,,t] %*% Fac[t,]
  }
  return(list(B=B, Fac=Fac, BF=BF))
}
