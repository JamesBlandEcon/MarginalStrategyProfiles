source("code/functions.R")
library(tidyverse)
library(viridis)
library(mnorm)
set.seed(42)

z<-halton(n=30) |> qnorm()


# Asymetric RPS

Payoffs<-rbind(c(1,0,4),c(9,1,0),c(0,16,1))
labels<-c("R","P","S")
print(paste(Payoffs,Payoffs |> t())) |> matrix(nrow=dim(Payoffs)[1])
nActions<-dim(Payoffs)[1]
rho0<-rep(log(1/nActions),nActions)


rSD<-c(0,0.2,0.4,0.8)

qre<-data.frame()

for (r in rSD) {

  X<-exp(log(0.5)+r*z)
  
  
  
  H<-function(x) {
    rho<-x[1:nActions]
    lambda<-x[1+nActions]
    
    p<-matrix(NA,nActions,length(X))
    
    for (xx in 1:length(X)) {
      p[,xx]<-softmax(lambda*(Payoffs^X[xx])%*%exp(rho))
    }
    
    h<-rep(NA,nActions)
    
    for (aa in 1:(nActions-1))  {
      h[aa]<-rho[aa+1]-rho[aa]-log(mean(p[aa+1,]))+log(mean(p[aa,]))
    }
    
    h[nActions]<-1-sum(exp(rho))
    
    h
  }
  
  H(c(rho0,1))
  
  Jac<-function(x) {
    rho<-x[1:nActions]
    lambda<-x[1+nActions]
    
    p<-exp(rho)
    
    Jac<-matrix(0,nActions,nActions+1)
    
    J_int1<-array(NA,dim=c(nActions,nActions,length(X)))
    J_int2<-array(NA,dim=c(nActions,1,length(X)))
    for (xx in 1:length(X)) {
      pay<-Payoffs^X[xx]
      u<-pay%*%p
      u_rho<-pay * (p %*% t(rep(1,nActions)))
      prob<-softmax(u)
      
      J_int1[,,xx]<-prob%*%t(rep(1,nActions))*(u_rho-(u_rho%*%prob)%*%t(rep(1,nActions)))
      
      J_int2[,1,xx]<-prob*(u-sum(prob*u))
    }
    
    mJ_int1<-apply(J_int1,FUN=mean,MARGIN=c(1,2))
    mJ_int2<-apply(J_int2,FUN=mean,MARGIN=c(1,2))
    for (aa in 1:(nActions-1)) {
      for (kk in 1:nActions) {
        Jac[aa,kk]<-(
          1*((aa+1)==kk)-1*(aa==kk)
          -lambda/exp(rho[aa+1])*mJ_int1[aa+1,kk]
          +lambda/exp(rho[aa])*mJ_int1[aa,kk]
        )
      }
      Jac[aa,nActions+1]<- (
        -1/exp(rho[aa+1])*mJ_int2[aa+1]
        +1/exp(rho[aa])*mJ_int2[aa]
      )
    }
    
    Jac[nActions,1:nActions]<--exp(rho)
    
    Jac
    
  }
  
  
  QRE<-PC(
    H,Jac,
    nActions,
    c(0,100),
    rho0,
    0.1,
    1e-4,
    1.1,
    1000,
    1e-4
  )
  
  colnames(QRE)<-c(paste0("rho",labels),"lambda")
  
  qre<-QRE |>
    data.frame() |>
    pivot_longer(cols = colnames(QRE)[1:nActions],
                 names_to="action",
                 values_to = "rho",
                 names_prefix = "rho"
    ) |>
    mutate(p=exp(rho),
           rSD=r) |>
    rbind(qre)

}
qre_wide<-qre |>
  pivot_wider(id_cols=c(lambda,rSD),names_from=action,values_from = p) |>
  mutate(rSDtext = paste0(rSD))

(
  ggplot(qre_wide,aes(x=R,y=P,group=rSD,linetype=rSDtext))
  +geom_path()
  +theme_bw()
  +coord_fixed()
  +labs(linetype="rSD")
)

list(qre=qre_wide,
        Payoffs = Payoffs
     ) |> 
  saveRDS("QREresults/SymmetricCRRA.rds")
