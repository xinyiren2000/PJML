rm(list = ls())
setwd("C:/Users/Administrator/Desktop/PJML/r")

source("gendata.r")
source("pjml.r")
source("others.r")
source("kernel.r")
source("numfac.r")
source("evaluation.r")

load_pkgs <- c("doParallel", "foreach")
sapply(load_pkgs, require, character = TRUE)

## returns names of functions loaded in current session
.export <- unclass(lsf.str())
.packages <- c("MASS", "mvtnorm", "base", "PDSCE")

sapply(.packages, require, character = TRUE)


## main function-------------------------------------------------------------------------
main_par <- function(HH=200, Ts, p, q, lam1.max=0.05, nods = 0,
                     alpha_f=0.3, s1=0.3, type_e = 'csd', sig=0.6, 
                     .export, .packages) {
  Sigma_e = gen_se(p, type_e = type_e, sig = sig)
  Ft = gen_f(Ts, q, alpha_f = alpha_f, f_0 = rep(0, q))
  B = gen_b(Ts, p, q, type_b = 'tv',s1=s1)
  BF = matrix(NA, nrow = Ts, ncol = p)
  for (t in 1:Ts) {
    BF[t,] = Ft[t,] %*% t(B[, , t])
  }
  if(nods==0)
  {
    lams.tv = CVPJML.tv(BF+gen_e(Ts, p, sigma_e = Sigma_e),q,lam1.max=lam1.max)
    write.csv(c(lam1=lams.tv$lam1,lam2=lams.tv$lam2),file = paste0("output/parameters/T",Ts,"p",p,"s",s1,"rhoe",sig,".csv"))
  } else{
    read.lams = read.csv(paste0("output/parameters/T",Ts,"p",p,"s",s1,"rhoe",sig,".csv"))
    lams.tv = list(lam1=read.lams[1,2],lam2=read.lams[2,2])
  }
  
  res <- foreach(
    h = 1:HH,
    .combine = rbind, .export = .export, .packages = .packages
  ) %dopar% {
    seed = 1234*(h+nods)
    set.seed(seed)
    E = gen_e(Ts, p, sigma_e = Sigma_e)
    Y = BF + E
    
    wls = local_pca(Y,q)
    apjml.tv = PJML.tv(Y=Y,q=q,lam1=lams.tv$lam1,lam2=lams.tv$lam2,ap=TRUE)
    pjml.tv = PJML.tv(Y=Y,q=q,lam1=lams.tv$lam1,lam2=lams.tv$lam2,ap=FALSE,maxit=3)
    
    ss <- c(seed =seed,
      wlsB=mDis(B,wls$B_tilde),pjmlB=mDis(B,pjml.tv$Bhat),apjmlB=mDis(B,apjml.tv$Bhat),
      wlsF=Dis(Ft,wls$pre_F),pjmlF=Dis(Ft,pjml.tv$Fhat),apjmlF=Dis(Ft,apjml.tv$Fhat)
    )
  }
  ## res is a matrix of HH rows and 24 columns.
  return(res)
}

#----sparse loading matrix----------------------------------------------
#Regular computation
TVCSD=NULL
for (Ts in c(50,100,200)) {
  for (p in c(50,100,200)) {
    for (sig in c(0,0.3,0.6)) {
      .export <- unclass(lsf.str())
      ## export packages if required
      .packages <- NULL
      
      library(foreach)
      
      PARTYPE <- "PC"
      
      if (PARTYPE == "PC") {
        ## cores used in local PC
        library(doParallel)
        cores <- detectCores() - 1
        cl <- makeCluster(cores)
        registerDoParallel(cl)
      } else if (PARTYPE == "HPC") {
        library(doMPI)
        cl <- startMPIcluster()
        registerDoMPI(cl)
        ## cleanup code, stop the cluster if an error occurs.
        options(error = function() {
          stopCluster(cl)
          options(error = NULL)
        })
      } else {
        stop("wrong parallel type!")
      }
      
      tvcsd = main_par(
        HH = 1000, Ts = Ts, p = p, q = 3, lam1.max=0.05, 
        alpha_f=0.3, s1=0.3, type_e = 'csd', sig = sig, 
        .export = .export, .packages = .packages
      )
      
      if (PARTYPE == "PC") {
        stopCluster(cl)
      } else if (PARTYPE == "HPC") {
        closeCluster(cl)
        mpi.quit()
      } else {
        stop("wrong parallel type!")
      }
      write.csv(tvcsd[,-1],file = paste0("output/s3/tvcsd_T", Ts,'_p',p,'_rhoe',sig, ".csv"))
      print(c(Ts,p,sig,colMeans(tvcsd[,-1])))
      TVCSD = rbind(TVCSD,c(Ts,p,sig,colMeans(tvcsd[,-1])),c(Ts,p,sig,apply(tvcsd[,-1],2,sd)))
    }
  }
}

write.csv(TVCSD,file = paste0("output/s3/tvcsd.csv"))


# Node-by-node computation
TVCSD=NULL
for (Ts in c(50,100,200)) {
  for (p in c(50,100,200)) {
    for (sig in c(0,0.3,0.6)) {
      simtime = 1000
      numnods = 10
      lnod = simtime/numnods
      tvcsd = NULL
      for(nods in ((1:numnods)-1)*lnod){
        print(nods/lnod)
        .export <- unclass(lsf.str())
        ## export packages if required
        .packages <- NULL
        
        library(foreach)
        
        PARTYPE <- "PC"
        
        if (PARTYPE == "PC") {
          ## cores used in local PC
          library(doParallel)
          cores <- detectCores() - 1
          cl <- makeCluster(cores)
          registerDoParallel(cl)
        } else if (PARTYPE == "HPC") {
          library(doMPI)
          cl <- startMPIcluster()
          registerDoMPI(cl)
          ## cleanup code, stop the cluster if an error occurs.
          options(error = function() {
            stopCluster(cl)
            options(error = NULL)
          })
        } else {
          stop("wrong parallel type!")
        }
        
        try = main_par(
          HH = lnod, Ts = Ts, p = p, q = 3, lam1.max=0.05, nods=nods,
          alpha_f=0.3, s1=0.3, type_e = 'csd', sig = sig, 
          .export = .export, .packages = .packages
        )
        
        if (PARTYPE == "PC") {
          stopCluster(cl)
        } else if (PARTYPE == "HPC") {
          closeCluster(cl)
          mpi.quit()
        } else {
          stop("wrong parallel type!")
        }
        print(c(Ts,p,sig,colMeans(try)))
        tvcsd = rbind(tvcsd,try)
        write.csv(tvcsd,file = paste0("output/s3/nods/tvcsd_T",Ts,'_p',p,'_rhoe',sig,"_nod",nods/lnod,".csv"))
      }
      write.csv(tvcsd[,-1],file = paste0("output/s3/tvcsd_T",Ts,'_p',p,'_rhoe',sig, ".csv"))
      print(c(Ts,p,sig,colMeans(tvcsd[,-1])))
      TVCSD = rbind(TVCSD,c(Ts,p,sig,colMeans(tvcsd[,-1])),c(Ts,p,sig,apply(tvcsd[,-1],2,sd)))
    }
  }
}

write.csv(TVCSD,file = paste0("output/s3/tvcsd.csv"))


#-----dense loading matrix--------------------------------------
# regular computation
TVCSD=NULL
for (Ts in c(50,100,200)) {
  for (p in c(50,100,200)) {
    for (sig in c(0,0.3,0.6,0.9)) {
      .export <- unclass(lsf.str())
      ## export packages if required
      .packages <- NULL
      
      library(foreach)
      
      PARTYPE <- "PC"
      
      if (PARTYPE == "PC") {
        ## cores used in local PC
        library(doParallel)
        cores <- detectCores() - 1
        cl <- makeCluster(cores)
        registerDoParallel(cl)
      } else if (PARTYPE == "HPC") {
        library(doMPI)
        cl <- startMPIcluster()
        registerDoMPI(cl)
        ## cleanup code, stop the cluster if an error occurs.
        options(error = function() {
          stopCluster(cl)
          options(error = NULL)
        })
      } else {
        stop("wrong parallel type!")
      }
      
      tvcsd = main_par(
        HH = 1000, Ts = Ts, p = p, q = 3, lam1.max=0.0, 
        alpha_f=0.3, s1=1, type_e = 'csd', sig = sig, 
        .export = .export, .packages = .packages
      )
      
      if (PARTYPE == "PC") {
        stopCluster(cl)
      } else if (PARTYPE == "HPC") {
        closeCluster(cl)
        mpi.quit()
      } else {
        stop("wrong parallel type!")
      }
      write.csv(tvcsd[,-1],file = paste0("output/dense/tvcsd_T",Ts,'_p',p,'_rhoe',sig, ".csv"))
      print(c(Ts,p,sig,colMeans(tvcsd[,-1])))
      TVCSD = rbind(TVCSD,c(Ts,p,sig,colMeans(tvcsd[,-1])),c(Ts,p,sig,apply(tvcsd[,-1],2,sd)))
    }
  }
}

write.csv(TVCSD,file = paste0("output/dense/tvcsd.csv"))

# Node-by-node computation for dense loading matrix is similar to that for sparse loading matrix.

