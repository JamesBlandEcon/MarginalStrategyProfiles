source("code/functions.R")
library(tidyverse)
library(mnorm)
set.seed(42)

nActions<-4

rho0<-rep(log(0.5),4)

x0<-c(rho0,1)

LSD<-c(0,1,2)


u1<-rbind(c(4,0),c(0,1))
u2<-rbind(c(0,1),c(1,0))

U<-rbind(cbind(matrix(0,2,2),u1),cbind(u2,matrix(0,2,2)))

storethis<-tibble()

for (lSD in LSD) {
  X<-exp(0+lSD*qnorm(halton(100)))
  
  
  H<-function(x) {
    rho<-x[1:nActions]
    lambda<-x[1+nActions]
    EU<-U%*% exp(rho)
    lEU<-EU %*% t(lambda*X)
    p<-matrix(NA,4,length(X))
    p[1:2,]<-exp(lEU[1:2,]-log(rep(1,2)%*%t(apply(exp(lEU[1:2,]),MARGIN=2,FUN=sum))))
    p[3:4,]<-exp(lEU[3:4,]-log(rep(1,2)%*%t(apply(exp(lEU[3:4,]),MARGIN=2,FUN=sum))))
    
    
    logmeanp<-log(apply(p,MARGIN=1,FUN=mean))
    h<-rep(NA,nActions)
    
    # Equilibrium constraints
    h[1]<-rho[2]-rho[1]-logmeanp[2]+logmeanp[1]
    h[2]<-rho[4]-rho[3]-logmeanp[4]+logmeanp[3]
    
    # Adding-up constraints
    h[3]<-1-sum(exp(rho[1:2]))
    h[4]<-1-sum(exp(rho[3:4]))
    
    
    h
  }
  
  jac<-function(x) {
    rho<-x[1:nActions]
    lambda<-x[1+nActions]
    EU<-U%*% exp(rho)
    lEU<-EU %*% t(lambda*X)
    p<-matrix(NA,4,length(X))
    p[1:2,]<-exp(lEU[1:2,]-log(rep(1,2)%*%t(apply(exp(lEU[1:2,]),MARGIN=2,FUN=sum))))
    p[3:4,]<-exp(lEU[3:4,]-log(rep(1,2)%*%t(apply(exp(lEU[3:4,]),MARGIN=2,FUN=sum))))
    prho<-matrix(NA,4,length(X))
    prho[1:2,]<-exp(lEU[1:2,]-rho[1:2]%*%t(rep(1,length(X)))-log(rep(1,2)%*%t(apply(exp(lEU[1:2,]),MARGIN=2,FUN=sum))))
    prho[3:4,]<-exp(lEU[3:4,]-rho[3:4]%*%t(rep(1,length(X)))-log(rep(1,2)%*%t(apply(exp(lEU[3:4,]),MARGIN=2,FUN=sum))))
    
    
    JAC<-matrix(NA,nActions,nActions+1)
    
    
    
    
    
    
    # part of the Jacobian corresponding to the derivative with respect to 
    # lambda
    eu1<-EU[1:2]%*%p[1:2,]
    eu2<-EU[3:4]%*%p[3:4,]
    hl<-apply(
      (prho*(EU%*%t(rep(1,length(X)))-rbind(eu1,eu1,eu2,eu1)))*rep(1,4)%*%t(X),
      MARGIN=1,FUN=mean
    )
    JAC[1,nActions+1]<-hl[2]-hl[1]
    JAC[2,nActions+1]<-hl[4]-hl[3]
    JAC[3:4,nActions+1]<-0
    
    
    
    # Part of the Jacobian corresponding to rho
    JAC[3,]<-c(-exp(rho[1:2]),0,0,0)
    JAC[4,]<-c(0,0,-exp(rho[3:4]),0)
    
    # Note here that we can separate out dUdrho from X, because we are 
    # multiplying through by X
    dUdrho<-U %*% diag(exp(rho))
    
    edudrho1<-t(dUdrho[1:2,])%*% p[1:2,]
    edudrho2<-t(dUdrho[3:4,])%*% p[3:4,]
    
    hrho1<-array(NA,c(2,4,length(X)))
    hrho2<-array(NA,c(2,4,length(X)))
    for (xx in 1:length(X)) {
      hrho1[,,xx]<-X[xx]*prho[1:2,xx] %*%t(rep(1,nActions))*(dUdrho[1:2,]-rep(1,2)%*%t(edudrho1[,xx]))
      hrho2[,,xx]<-X[xx]*prho[3:4,xx]%*%t(rep(1,nActions))*(dUdrho[3:4,]-rep(1,2)%*%t(edudrho2[,xx]))
    }
    
    HRHO1<-apply(hrho1,MARGIN=c(1,2),FUN=mean)
    HRHO2<-apply(hrho2,MARGIN=c(1,2),FUN=mean)
    
    JAC[1,1:nActions]<-lambda*(-HRHO1[2,]+HRHO1[1,])+c(-1,1,0,0)
    JAC[2,1:nActions]<-lambda*(-HRHO2[2,]+HRHO2[1,])+c(0,0,-1,1)
    
    JAC
  }
  
  QRE<-PC(
    H,jac,
    nActions,
    c(0,100),
    rho0,
    0.001,
    1e-4,
    1.001,
    1000,
    1e-4
  )
  
  colnames(QRE)<-c("rhoU","rhoD","rhoL","rhoR","lambda")
  
  qre<-QRE |> 
    data.frame() |>
    mutate(
      U = exp(rhoU),
      D = exp(rhoD),
      L = exp(rhoL),
      R = exp(rhoR ),
      lSD = lSD
    )
  
  storethis<-rbind(storethis,qre)
}

(
  ggplot(storethis,aes(x=U,y=L,linetype=as.factor(lSD)))
  
  +geom_path()
  +theme_bw()
  #+coord_fixed()
  +coord_fixed()
)

storethis |> saveRDS("QREresults/Ochs1995.rds")
