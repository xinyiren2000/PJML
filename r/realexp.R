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
.packages <- c("MASS", "mvtnorm", "base", "PDSCE","fbi")

sapply(.packages, require, character = TRUE)


# install.packages("stats")
# install.packages("readr")
# install.packages("pracma")
# devtools::install_github("cykbennie/fbi")
fbi.usme=fbi::fredmd(file="us_macroeconometrics/2024-12-m.csv")
#sum(is.na(fbi.usme[399:735,]))
#Ydat = scale(as.matrix(fbi.usme[399:735,-1])) #1992/3-2020/3
Ydat = scale(as.matrix(fbi.usme[399:732,-1])) #1992/3-2019/12


roll_par <- function(Y, q=NULL, width=100, step=1, lam1.max=0.05) {
  Ts = nrow(Y)
  p = ncol(Y)
  HH = Ts-width-step+1
  
  if(is.null(q)){
    q1 = select_q(Y[1:width,])
    lams.tv = CVPJML.tv(Y[1:width,],q1,lam1.max=lam1.max)
  }else{
    lams.tv = CVPJML.tv(Y[1:width,],q,lam1.max=lam1.max)
  }
  
  res <- foreach(
    h = 1:HH,
    .combine = rbind, .export = .export, .packages = .packages
  ) %dopar% {
    #pre.name = c("INDPRO","RPI","HWI","UNRATE","DPCERA3M086SBEA","ISRATIOx","WPSFD49207","CUSR0000SA0L2")
    Ytrain = Y[h:(h+width-1),]
    Ytest = Y[h+width+step-1,,drop=FALSE]
    if(is.null(q)){
      q = select_q(Ytrain)
    }
    
    wls = local_pca(Ytrain,q)
    apjml.tv = PJML.tv(Y=Ytrain,q=q,lam1=lams.tv$lam1,lam2=lams.tv$lam2,ap=TRUE)
    pjml.tv = PJML.tv(Y=Ytrain,q=q,lam1=lams.tv$lam1,lam2=lams.tv$lam2,ap=FALSE,maxit=3)
    
    Awls = t(Ytrain[2:width,])%*%wls$pre_F[1:(width-1),]%*%solve(t(wls$pre_F[1:(width-1),])%*%wls$pre_F[1:(width-1),])
    Ywls = wls$pre_F[width,,drop=F] %*% t(Awls)
    
    Apjml = t(Y[2:width,])%*%pjml.tv$Fhat[1:(width-1),]%*%solve(t(pjml.tv$Fhat[1:(width-1),])%*%pjml.tv$Fhat[1:(width-1),])
    Ypjml = pjml.tv$Fhat[width,,drop=F] %*% t(Apjml)
    
    Aapjml = t(Y[2:width,])%*%apjml.tv$Fhat[1:(width-1),]%*%solve(t(apjml.tv$Fhat[1:(width-1),])%*%apjml.tv$Fhat[1:(width-1),])
    Yapjml = apjml.tv$Fhat[width,,drop=F] %*% t(Aapjml)
    
    
    ss <- c(
      q=q,
      ((Ytest-Ywls)^2), ((Ytest-Ypjml)^2), ((Ytest-Yapjml)^2)
    )
  }
  return(res)
}



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

try = roll_par(Y=Ydat[,],q=NULL,width=100)

if (PARTYPE == "PC") {
  stopCluster(cl)
} else if (PARTYPE == "HPC") {
  closeCluster(cl)
  mpi.quit()
} else {
  stop("wrong parallel type!")
}
results = t(matrix(colMeans(try[,-1]),ncol = 3))
colnames(results) = colnames(Ydat)
rownames(results) = c('wls','pjml','apjml')
better.results = results[,results[1,]>results[2,]]
print(better.results)

p = ncol(Ydat)
write.csv(try[,1:(p+1)],file = paste0("real data/qroll/wls.csv"))
write.csv(try[,c(1,(p+2):(2*p+1))],file = paste0("real data/qroll/pjml.csv"))
write.csv(try[,c(1,(2*p+2):(3*p+1))],file = paste0("real data/qroll/apjml.csv"))
write.csv(results,file = paste0("real data/qroll/results.csv"))
write.csv(better.results,file = paste0("real data/qroll/better.csv"))



for(q in 2:5){
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
  
  try = roll_par(Y=Ydat[,],q=q,width=100,step=1,lam1.max = 0.05)
  
  if (PARTYPE == "PC") {
    stopCluster(cl)
  } else if (PARTYPE == "HPC") {
    closeCluster(cl)
    mpi.quit()
  } else {
    stop("wrong parallel type!")
  }
  results = t(matrix(colMeans(try[,-1]),ncol = 3))
  colnames(results) = colnames(Ydat)
  rownames(results) = c('wls','pjml','apjml','pca','wpll')
  better.results = results[,results[1,]>results[2,]]
  print(better.results)
  
  p = ncol(Ydat)
  write.csv(try[,2:(p+1)],file = paste0("real data/q",q,"/wls.csv"))
  write.csv(try[,(p+2):(2*p+1)],file = paste0("real data/q",q,"/pjml.csv"))
  write.csv(try[,(2*p+2):(3*p+1)],file = paste0("real data/q",q,"/apjml.csv"))
  write.csv(results,file = paste0("real data/q",q,"/results.csv"))
  write.csv(better.results,file = paste0("real data/q",q,"/better.csv"))
}




