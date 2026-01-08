select_q <- function(Y, q_max=10, eps = 0.01, method = "BIC"){
  p = ncol(Y)
  Ts = nrow(Y)
  h = (2.35/sqrt(12)) * Ts^(-1/5) * p^(-1/10)
  Ys = array(data=NA, dim=c(Ts,p,Ts))
  if(method=='ratio'){
    Q = rep(NA,Ts)
    for(s in 1:Ts){ 
      for (t in 1:Ts){
        Ys[t,,s] = Y[t,]*k(Ts,p,t,s)^(1/2)
      }
      # kernel-weighted sample covariance matrix
      lambdas = eigen(Ts^(-1) * t(Ys[,,s])%*%Ys[,,s])$values
      rlambda = rep(NA,q_max)
      for(q in 1:q_max){
        if(abs(lambdas[q]/(p)) <= eps){
          rlambda[q] = 1
        }else{
          rlambda[q] = lambdas[q+1]/lambdas[q]
        }
      }
      Q[s] = which.min(rlambda)
    }
    qopt = max(Q) #max(Q[round((1:9)*(Ts/10))])
  }else if(method=='BIC'){
    Q = rep(NA,q_max)
    for (q in 2:q_max) {
      BF = matrix(data=NA, nrow=Ts, ncol=p)
      for (s in 1:Ts){
        for (t in 1:Ts){
          Ys[t,,s] = Y[t,]*k(Ts,p,t,s)^(1/2)
        }
        Fs_tilde = sqrt(Ts)*as.matrix(eigen(Ys[,,s]%*%t(Ys[,,s]))$vectors)[,1:q]
        Bs_tilde = Ts^(-1) * t(Ys[,,s])%*%Fs_tilde
        Bs_check = (p*Ts)^(-1) * t(Ys[,,s])%*%Ys[,,s]%*%Bs_tilde
        pre_Fscheck = Y[s,]%*%Bs_check %*%solve(t(Bs_check)%*%Bs_check)  
        BF[s,] = pre_Fscheck%*%t(Bs_check)
      }
      Q[q] = log(mean((Y-BF)^2)) + q*((p+Ts*h)/(p*Ts*h))*log((p*Ts*h)/(p+Ts*h))
    }
    qopt = which.min(Q[-1])+1
    
  }
  return(qopt)
}
