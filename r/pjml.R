TSPCA <- function(Y, q, lam1=NULL, lam2=NULL, 
                  cov.tol=1e-2, cov.maxit=100){
  p = ncol(Y)
  Ts = nrow(Y)
  pdiag = (Ts <= p)
  
  old.F = sqrt(Ts)*as.matrix(eigen(Y%*%t(Y))$vectors[,1:q])
  old.B=t(Y)%*%old.F / Ts
  BF = old.F%*%t(old.B)
  res = Y-BF
  residual.cov = crossprod(res) / Ts
  siom=PDSCE::pdsoft(s=residual.cov,lam=lam2,init='soft',tolin = cov.tol, 
                     tolout = cov.tol, maxitin = cov.maxit, maxitout = cov.maxit)
  Sigma=siom$sigma
  om=siom$omega
  
  Z = Ts^(-1/2) * Y%*%om%*%old.B
  W = as.matrix(eigen(t(Z)%*%Z)$vectors[,1:q])
  V = as.matrix(eigen(Z%*%t(Z))$vectors[,1:q])
  new.F = sqrt(Ts) * V%*%t(W)
  for(j in 1:q){
    if(sign(old.F[1,j])!=sign(new.F[1,j]))
      new.F[,j] = new.F[,j] * (-1)
  }
  
  refit=MRCE::mrce(Y=Y, X=new.F, lam2=lam1, method="fixed.omega", omega=om) 
  new.B0 = t(refit$Bhat)
  svd.B = svd(new.B0)
  new.B = new.B0%*%svd.B$v
  
  new.F0 = new.F
  new.F = new.F0%*%svd.B$v
  
  BF = new.F%*%t(new.B)
  
  return(list(Bhat=new.B,Fhat=new.F,BF=BF,sigma=Sigma,omega=om))
}

#tspca = TSPCA(Y,q,lam1=0.05,lam2=0.1)

TSPCA.tv <- function(Y, q, lam1=NULL, lam2=NULL, 
                     cov.tol=1e-2, cov.maxit=100){
  p = ncol(Y)
  Ts = nrow(Y)
  h = (2.35/sqrt(12)) * Ts^(-1/5) * p^(-1/10)
  Ys = array(data=NA, dim=c(Ts,p,Ts))
  Bhat = array(data=NA, dim=c(p,q,Ts))
  Fhat = matrix(data=0, nrow=Ts, ncol=q)
  BF = matrix(data=NA, nrow=Ts, ncol=p)
  
  for (s in 1:Ts){
    for (t in 1:Ts){
      Ys[t,,s] = Y[t,]*k(Ts,p,t,s)^(1/2)
    }
    fits = TSPCA(Ys[,,s],q,lam1=lam1,lam2=lam2,cov.tol=cov.tol, cov.maxit=cov.maxit)
    Bhat[,,s] = fits$Bhat
    if(s>1){
      for(j in 1:q){
        if(sign(Bhat[1,j,1])!=sign(Bhat[1,j,s]))
          Bhat[,j,s] = Bhat[,j,s] * (-1)
      }
    }
    # pre_F[s,] <- tryCatch(
    #   solve(t(B_hat[,,s])%*%B_hat[,,s]) %*% t(B_hat[,,s])%*%Y[s,],
    #   error = function(e) solve(t(B_hat[,,s])%*%B_hat[,,s]+lam1*diag(1,q)) %*% t(B_hat[,,s])%*%Y[s,]
    # )
    Fs = fits$Fhat
    Fs1 = matrix(NA,Ts,q)
    for(t in 1:Ts){
      Fs1[t,] = Fs[t,] * k(Ts,p,t,s)^(1/2)
    }
    Fhat = Fhat+Fs1*Ts^(-1) 
    # BF[s,] = B_hat[,,s]%*%pre_F[s,]
  }
  for(s in 1:Ts){
    BF[s,] = Bhat[,,s]%*%Fhat[s,]
  }
  return(list(Bhat=Bhat, Fhat=Fhat, BF=BF))
}


#tspca.tv = TSPCA.tv(Y,q,lam1=0.05,lam2=0.1)




joint_PCA <- function(Y, q, lam1=NULL, lam2=NULL, penal='mcp', diagBB = FALSE,
                      cov.tol=1e-2, cov.maxit=100, maxit=100, eps=1e-2){
  aw_mcp <- function(u, lambda, a=2.5) {
    if (a <= 1) {
      stop("a should be larger than 1.")
    }
    pen_prime_abs <- pmax((a * lambda - abs(u)) / a, 0)
    return(pen_prime_abs)
  }
  
  p = ncol(Y)
  Ts = nrow(Y)
  pdiag = (Ts <= p)
  
  old.F = sqrt(Ts)*as.matrix(eigen(Y%*%t(Y))$vectors[,1:q])
  old.B=t(Y)%*%old.F / Ts
  BF = old.F%*%t(old.B)
  res = Y-BF
  residual.cov = crossprod(res) / Ts
  sigma=diag(diag(residual.cov))
  om=diag(1/diag(residual.cov))
  
  iterating = TRUE
  k = 0
  
  siom=PDSCE::pdsoft(s=residual.cov,lam=lam2,init='soft',tolin = cov.tol, 
                     tolout = cov.tol, maxitin = cov.maxit, maxitout = cov.maxit)
  Sigma=siom$sigma
  om=siom$omega 
  
  Z = Ts^(-1/2) * Y%*%om%*%old.B
  W = as.matrix(eigen(t(Z)%*%Z)$vectors[,1:q])
  V = as.matrix(eigen(Z%*%t(Z))$vectors[,1:q])
  new.F = sqrt(Ts) * V%*%t(W)
  for(j in 1:q){
    if(sign(old.F[1,j])!=sign(new.F[1,j]))
      new.F[,j] = new.F[,j] * (-1)
  }
  #fdist = base::norm(old.F-new.F,type='F')^2/(Ts*q)
  old.F = new.F
  
  while (iterating) {
    k=k+1
    
    ftf=crossprod(old.F)
    omresf = om%*%t(res)%*%old.F 
    new.B = matrix(NA,p,q)
    for (i in 1:p) {
      for(j in 1:q){
        mainpart = old.B[i,j]+omresf[i,j]/(om[i,i]*ftf[j,j])
        if(penal=='mcp'){
          pen = aw_mcp(u=old.B[i,j], lambda=lam1)
        }else if(penal=='l1'){
          pen = lam1
        }
        new.B[i,j] = sign(mainpart)*max(abs(mainpart)-Ts*pen/(om[i,i]*ftf[j,j]),0)
      }
    }
    bdist = base::norm(old.B-new.B,type='F')^2/(p*q)
    old.B = new.B
    
    iterating = (bdist>eps) & (k <= maxit)
    BF = old.F%*%t(old.B)
    res = Y-BF
    #residual.cov = crossprod(res) / Ts
  }
  
  return(list(Bhat=old.B,Fhat=old.F,BF=BF,omega=om,num.it=k))
}


PJML.bf <- function(Y,q,si=NULL,lam1,lam2=NULL,gamma=1.5,
                    maxit=100,eps=1e-08, maxit.v=100,eps.v=1e-08,
                    cov.tolin = 1e-08, cov.tolout = 1e-08, 
                    cov.maxitin = 10000, cov.maxitout = 1000){
  Ts = nrow(Y)
  p = ncol(Y)
  svd.y = svd(Y)
  
  U.old = svd.y$u[,1:q,drop=FALSE]
  D.old = diag(svd.y$d[1:q])
  V.old = svd.y$v[,1:q,drop=FALSE]
  
  mu.old = 1e-03
  G.old = matrix(0,nrow=p,ncol=q)
  
  B.old = sqrt(1/Ts) * V.old%*%D.old
  F.old = sqrt(Ts) * U.old
  BF = F.old%*%t(B.old)
  res = Y-BF
  res.cov = crossprod(res) / Ts
  if(is.null(si)){
    if(is.null(lam2)){
      siom=PDSCE::pdsoft.cv(x=res,init="soft",tolin=cov.tolin,tolout=cov.tolout,
                            maxitin=cov.maxitin,maxitout=cov.maxitout)
      lam2 = siom$best.lam
    }else{
      siom=PDSCE::pdsoft(s=res.cov,lam=lam2,init='soft',
                         tolin=cov.tolin,tolout=cov.tolout,
                         maxitin=cov.maxitin,maxitout=cov.maxitout)
    }
    si=siom$sigma
    om=siom$omega
  }else{
    om=solve(si)
  }
  
  iterating = TRUE
  k = 0
  loss = NULL
  
  while (iterating) {
    k=k+1
    
    # U-step
    Z = Y%*%om%*%V.old%*%D.old
    svd.z =svd(Z)
    U.new = svd.z$u %*% t(svd.z$v)
    
    # V-step
    V.new = vstep.wop(om=om,Y=Y,U=U.new,D=D.old,V0=V.old,B=B.old,G=G.old,mu=mu.old,
                      maxit=maxit.v,eps=eps.v)$V
    
    # D-step
    C = sqrt(Ts)*(B.old-(1/mu.old)*G.old)
    D.new = diag(NA,q)
    for (j in 1:q) {
      uj = U.new[,j,drop=FALSE]
      vj = V.new[,j,drop=FALSE]
      cj = C[,j,drop=FALSE]
      dm = t(uj)%*%Y%*%om%*%vj+mu.old*t(vj)%*%cj
      dd = t(vj)%*%om%*%vj+mu.old
      D.new[j,j] = dm/dd
    }
    
    # B-step 
    W = sqrt(1/Ts)*V.new%*%D.new + (1/mu.old)*G.old
    B.new = matrix(NA,nrow=p,ncol=q)
    for(i in 1:p){
      for(j in 1:q){
        wij = W[i,j]
        B.new[i,j] = sign(wij)*max(abs(wij)-lam1/(2*mu.old),0)
      }
    }
    #F.new = sqrt(Ts) * U.new
    
    #G-step 
    G.new = G.old + mu.old*(sqrt(1/Ts)*V.new%*%D.new-B.new)
    mu.new = gamma*mu.old
    
    bdist = base::norm(B.old-B.new,type='F')^2
    udist = base::norm(U.old-U.new,type='F')^2
    #fdist = base::norm(F.old-F.new,type='F')^2/(Ts*q)
    loss.k = bdist+udist
    loss = rbind(loss,c(k=k,sumdis=loss.k,bdis=bdist,udis=udist))
    
    iterating = (loss.k>=eps) & (k < maxit)
    B.old = B.new
    U.old = U.new
    V.old = V.new
    D.old = D.new
    G.old = G.new
    mu.old = mu.new
  }
  return(list(Bhat=B.old,Fhat=sqrt(Ts)*U.old,N.it=k,loss=loss))
}

vstep.wop <- function(om,Y,U,D,V0,B,G,mu,maxit=100,eps=0.001){
  Ts = nrow(Y)
  p = ncol(Y)
  rho2 = max(eigen(om)$values)
  XX = rho2*diag(1,p)-om
  V.old = V0
  
  iterating = TRUE
  m = 0
  v.dis = NULL
  
  while (iterating) {
    m=m+1
    
    A = (om%*%t(Y)%*%U + sqrt(Ts)*(mu*B-G) + XX%*%V.old%*%D) %*% D
    svd.a = svd(A)
    V.new = svd.a$u %*% t(svd.a$v)
    
    dist = base::norm(V.old-V.new,type='F')^2
    v.dis = rbind(v.dis,c(m=m,dist=dist))
    V.old = V.new
    
    iterating = (dist>=eps) & (m < maxit)
  }
  return(list(V=V.old, N.it=m, v.dis=v.dis))
}

PJML <- function(Y,q,lam1,lam2=NULL,gamma=1.5,ap=FALSE,maxit=100,eps=1e-02,
                 maxit.in=100,eps.in=1e-08, maxit.v=100,eps.v=1e-08,
                 cov.tolin = 1e-08, cov.tolout = 1e-08, 
                 cov.maxitin = 10000, cov.maxitout = 1000){
  Ts = nrow(Y)
  p = ncol(Y)
  
  fit = PJML.bf(Y=Y,q=q,si=NULL,lam1=lam1,lam2=lam2,gamma=gamma,
                maxit=maxit.in,eps=eps.in, maxit.v=maxit.v,eps.v=eps.v,
                cov.tolin = cov.tolin, cov.tolout = cov.tolout, 
                cov.maxitin = cov.maxitin, cov.maxitout = cov.maxitout)
  B.old = fit$Bhat
  F.old = fit$Fhat
  kk = fit$N.it
  loss = fit$loss
  
  if(ap==FALSE){
    BF = F.old%*%t(B.old)
    res = Y-BF
    res.cov = crossprod(res) / Ts
    if(is.null(lam2)){
      siom=PDSCE::pdsoft.cv(x=res,init="soft",tolin=cov.tolin,tolout=cov.tolout,
                            maxitin=cov.maxitin,maxitout=cov.maxitout)
      lam2 = siom$best.lam
    }else{
      siom=PDSCE::pdsoft(s=res.cov,lam=lam2,init='soft',
                         tolin=cov.tolin,tolout=cov.tolout,
                         maxitin=cov.maxitin,maxitout=cov.maxitout)
    }
    si=siom$sigma
    om=siom$omega
    
    iterating = TRUE
    kk = 0
    loss = NULL
    
    while (iterating) {
      kk=kk+1
      
      fit = PJML.bf(Y=Y,q=q,si=si,lam1=lam1,lam2=lam2,gamma=gamma,
                    maxit=maxit.in,eps=eps.in, maxit.v=maxit.v,eps.v=eps.v)
      B.new = fit$Bhat
      F.new = fit$Fhat
      
      bdist = base::norm(B.old-B.new,type='F')^2
      fdist = base::norm(F.old-F.new,type='F')^2/Ts
      #fdist = base::norm(F.old-F.new,type='F')^2/(Ts*q)
      loss.k = bdist+fdist
      loss = rbind(loss,c(kk=kk,sumdis=loss.k,bdis=bdist,fdis=fdist,nit=fit$N.it))
      
      iterating = (loss.k>=eps) & (kk < maxit)
      B.old = B.new
      F.old = F.new
      
      BF = F.old%*%t(B.old)
      res = Y-BF
      res.cov = crossprod(res) / Ts
      siom=PDSCE::pdsoft(s=res.cov,lam=lam2,init='soft',
                         tolin=cov.tolin,tolout=cov.tolout,
                         maxitin=cov.maxitin,maxitout=cov.maxitout)
      si=siom$sigma
      om=siom$omega
    }
  }
  return(list(Bhat=B.old,Fhat=F.old,N.it=kk,loss=loss))
}

CVPJML <- function(Y,q,gamma=1.5,kfold=5,
                   lam1.max=NULL,lam2=NULL,nlam=5,lam.min.ratio= 0.01,
                   maxit=100,eps=1e-08, maxit.v=100,eps.v=1e-08,
                   cov.tolin = 1e-08, cov.tolout = 1e-08, 
                   cov.maxitin = 10000, cov.maxitout = 1000){
  Ts = nrow(Y)
  p = ncol(Y)
  svd.y = svd(Y)
  U.old = svd.y$u[,1:q,drop=FALSE]
  D.old = diag(svd.y$d[1:q])
  V.old = svd.y$v[,1:q,drop=FALSE]
  B.old = sqrt(1/Ts) * V.old%*%D.old
  F.old = sqrt(Ts) * U.old
  BF = F.old%*%t(B.old)
  res = Y-BF
  res.cov = crossprod(res) / Ts
  if(is.null(lam2)){
    siom=PDSCE::pdsoft.cv(x=res,init="soft",tolin=cov.tolin,tolout=cov.tolout,
                          maxitin=cov.maxitin,maxitout=cov.maxitout)
    lam2 = siom$best.lam
  }else{
    siom=PDSCE::pdsoft(s=res.cov,lam=lam2,init='soft',
                       tolin=cov.tolin,tolout=cov.tolout,
                       maxitin=cov.maxitin,maxitout=cov.maxitout)
  }
  si=siom$sigma
  om=siom$omega
  
  if(lam1.max<=0){
    lam1=0
  }else{
    if(is.null(lam1.max)){
      lam1.max = sqrt(log(Ts)/Ts) 
    }
    lam1.min = lam.min.ratio * lam1.max
    # calculate grid of lambda values
    lam1.vec = 10^seq(log10(lam1.min), log10(lam1.max), length = nlam+1)[2:nlam+1]
    
    err.mat = matrix(NA,kfold,length(lam1.vec))
    ind = 1:Ts
    for(k in 1:kfold){
      foldind = ind[ (1+floor((k-1)*Ts/kfold)):floor(k*Ts/kfold) ]
      Y.tr=Y[-foldind,]
      Y.va=Y[foldind,,drop=FALSE ]
      #F.va=pca$pre_F[foldind,,drop=FALSE ]
      #F.va = sqrt(nrow(Y.va))*as.matrix(eigen(Y.va%*%t(Y.va))$vectors[,1:q])
      for(i in 1:length(lam1.vec)){
        fit.tr = PJML.bf(Y=Y.tr,q=q,si=si,lam1=lam1.vec[i],lam2=lam2,gamma=gamma,
                         maxit=maxit,eps=eps, maxit.v=maxit.v,eps.v=eps.v)
        B.tr = fit.tr$Bhat
        F.tr = fit.tr$Fhat
        fit.va = PJML.bf(Y=Y.va,q=q,si=si,lam1=lam1.vec[i],lam2=lam2,gamma=gamma,
                         maxit=maxit,eps=eps, maxit.v=maxit.v,eps.v=eps.v)
        B.va = fit.va$Bhat
        F.va = fit.va$Fhat
        err.mat[k,i] = norm(F.va%*%t(B.tr)-Y.va,type = 'F')^2+norm(F.tr%*%t(B.va)-Y.tr,type='F')^2
      }
    }
    ## Eliminate saturated lambda values, if any
    el_ind <- which(apply(is.finite(err.mat), 2, all))
    err.mat <- err.mat[, el_ind, drop=FALSE]
    lam1.vec <- lam1.vec[el_ind]
    cve = apply(err.mat,2,mean)
    lam1 = lam1.vec[which.min(cve)]
  }
  return(list(lam1=lam1,lam2=lam2))
}

PJML.tv <- function(Y,q,lam1,lam2=NULL,gamma=1.5,
                    ap=FALSE,bm=TRUE,maxit=100,eps=1e-02,
                    maxit.in=100,eps.in=1e-08, maxit.v=100,eps.v=1e-08,
                    cov.tolin = 1e-08, cov.tolout = 1e-08, 
                    cov.maxitin = 10000, cov.maxitout = 1000){
  p = ncol(Y)
  Ts = nrow(Y)
  h = (2.35/sqrt(12)) * Ts^(-1/5) * p^(-1/10)
  Ys = array(data=NA, dim=c(Ts,p,Ts))
  B = array(data=NA, dim=c(p,q,Ts))
  preF = matrix(data=0, nrow=Ts, ncol=q)
  BF = matrix(data=NA, nrow=Ts, ncol=p)
  N.it = rep(NA,Ts)
  
  for (s in 1:Ts){
    for (t in 1:Ts){
      Ys[t,,s] = Y[t,]*k(Ts=Ts,p=p,t=t,s=s,boundary.modify=bm)^(1/2)
    }
    fit = PJML(Y=Ys[,,s],q=q,lam1=lam1,lam2=lam2,gamma=gamma,ap=ap,maxit=maxit,eps=eps,
               maxit.in=maxit.in,eps.in=eps.in,maxit.v=maxit.v,eps.v=eps.v,
               cov.tolin=cov.tolin,cov.tolout=cov.tolout,
               cov.maxitin=cov.maxitin,cov.maxitout=cov.maxitout)
    B[,,s] = fit$Bhat
    # if(s>1){
    #   for(j in 1:q){
    #     if(sign(B[1,j,1])!=sign(B[1,j,s]))
    #       B[,j,s] = B[,j,s] * (-1)
    #   }
    # }
    # pre_F[s,] <- tryCatch(
    #   solve(t(B_hat[,,s])%*%B_hat[,,s]) %*% t(B_hat[,,s])%*%Y[s,],
    #   error = function(e) solve(t(B_hat[,,s])%*%B_hat[,,s]+lam1*diag(1,q)) %*% t(B_hat[,,s])%*%Y[s,]
    # )
    Fs = fit$Fhat
    Fs1 = matrix(NA,Ts,q)
    for(t in 1:Ts){
      Fs1[t,] = Fs[t,] * k(Ts=Ts,p=p,t=t,s=s,boundary.modify=bm)^(1/2)
    }
    preF = preF+Fs1*Ts^(-1) 
    # BF[s,] = B_hat[,,s]%*%pre_F[s,]
    N.it[s] = fit$N.it
  }
  for(s in 1:Ts){
    BF[s,] = B[,,s]%*%preF[s,]
  }
  return(list(Bhat=B, Fhat=preF, BF=BF, num.it=N.it))
}

CVPJML.tv <- function(Y,q,gamma=1.5,bm=TRUE,kfold=5,
                      lam1.max=NULL,lam2=NULL,nlam=5,lam.min.ratio= 0.01,
                      maxit=100,eps=1e-08, maxit.v=100,eps.v=1e-08,
                      cov.tolin = 1e-08, cov.tolout = 1e-08, 
                      cov.maxitin = 10000, cov.maxitout = 1000){
  p = ncol(Y)
  Ts = nrow(Y)
  best.lam1=0
  best.lam2=0
  for (s in round((1:9)*(Ts/10))){
    Ys = matrix(NA,Ts,p)
    for (t in 1:Ts){
      Ys[t,] = Y[t,,drop=F]*k(Ts=Ts,p=p,t=t,s=s,boundary.modify=bm)^(1/2)
    }
    lams = CVPJML(Y=Ys,q=q,kfold=kfold,lam1.max=lam1.max,lam2=lam2, 
                  nlam=nlam,lam.min.ratio=lam.min.ratio, 
                  maxit.v=maxit.v,eps.v=eps.v,
                  cov.tolin=cov.tolin,cov.tolout=cov.tolout,
                  cov.maxitin=cov.maxitin,cov.maxitout=cov.maxitout)
    #print(lams)
    best.lam1=best.lam1+lams$lam1/9
    best.lam2=best.lam2+lams$lam2/9
  }  
  return(list(lam1=best.lam1,lam2=best.lam2))
}
