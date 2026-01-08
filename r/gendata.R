gen_f <- function(Ts, q=2, alpha_f=0.3, f_0=c(0,0)){
  set.seed(12345)
  Ft = matrix(data=NA, nrow=q, ncol=Ts+1)
  Ft[,1] = f_0
  for (t in 1:Ts){
    for (j in 1:q){
      Ft[j,t+1] = alpha_f*Ft[j,t] + rnorm(1,mean=0,sd=sqrt(1-alpha_f^2))
    }
  }
  return(t(Ft[,-1]))
}
# gen_f(20)

#library(MASS)
#library(mvtnorm)
gen_b <- function(Ts, p, q=2, type_b='tv', s1=0.5){
  set.seed(123456)
  
  if (type_b == "static"){
    B = matrix(rnorm(p*q), p, q)*matrix(rbinom(p*q,1,s1), p, q)
  } else if (type_b == "tv"){
    B = array(data=NA, dim=c(p,q,Ts))
    for (i in 1:p){
      for (j in 1:q){
        bij  =  rnorm(1,mean=0,sd=1)*rbinom(1,1,s1)  #runif(1,-1,1)
        for (t in 1:Ts)
          #B[i,j,t] = runif(1,-2,2) + runif(1,-4,4)/(1+exp(-2*(10*t/Ts-runif(1,2,7))))
          B[i,j,t] = bij/(1+exp(-2*(10*t/Ts-(5*i/p+2)))) #bij + 0.1*exp(-0.7+3.5*(t/Ts)^(j-1)) 
        # 2*(t/Ts)^(j-1) + exp(-16*(t/Ts-0.5)^2) - 1 
        #* (t/Ts)^(j-1)
      }
    }
  }
  return(B)
}


# gen_b(20, 10, q=2, type_b='tv')



gen_se <- function(p,type_e='iid',sig=0.5){
  if (type_e=='csd'){
    # cross sectional dependence
    Sigma_e = matrix(NA, nrow=p, ncol=p)
    for (i in 1:p){
      for (j in 1:p){
        Sigma_e[i,j]=sig^(abs(i-j))
      }
    }
  }else if(type_e=='csh'){
    # cross sectional heteroskedasticity
    Sigma_e = diag(runif(p,1-sig,1+sig)^2)
  }else{
    # iid
    Sigma_e = diag(1,p)
  }
  return(Sigma_e)
}


gen_e <- function(Ts, p, sigma_e){
  E = MASS::mvrnorm(Ts, mu=rep(0,p), Sigma=sigma_e)
  return(E)
}
# gen_e(20,10,gen_se(10,'csd'))



#Y=MASS::mvrnorm(20,rep(0,10),diag(1,10))

