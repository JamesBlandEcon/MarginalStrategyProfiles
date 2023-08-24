# Needed for the Moore-Penrose inverse
library(MASS)

PC<-function(
  H, jac ,
  # Total number of actions in the game
  nActions,
  # range of lambda to compute QRE for. First element is the starting point. Second element can be above or below this
  lambda_span, # c(lambda0,lambda_end)
  # Starting poring for log probabilities
  x0, # rho0
  # first step size
  first_step,
  # Minimum step size
  min_step,
  # Maximum deceleration of step size
  max_decel,
  # Maimum number of iterations for corrector step
  maxiter,
  # Tolerance for corrector step (a step size)
  tol
) {
  
  
  # Maximum distance to curve
  max_dist <- 0.4
  # Maximum contraction rate in corrector
  max_contr <- 0.6
  # Perturbation to avoid cancellation in calculating contraction rate
  eta <- 0.1
  
  # Initial conditions
  x<-c(x0,lambda_span[1])
  QRE<-rbind(c(),x)
  
  
  
  # The last column of Q in the QR decomposition is the tangent vector
  t <- (qr(jac(x) |> t())  |> qr.Q(complete=TRUE))[,nActions+1]
  
  h<-first_step
  
  # Set orientation of curve, so we are tracing into the interior of `t_span`
  omega<-sign(t[nActions+1])*sign(lambda_span[2]-lambda_span[1])
  
  
  while ((x[nActions+1]>=min(lambda_span)) & (x[nActions+1]<=max(lambda_span))) {
    accept <-TRUE
    
    if (abs(h)<= min_step) {
      # Stepsize below minimum tolerance; terminate.
      success <-FALSE
      print(paste("Stepsize",abs(h),"less than minimum",min_step))
      break
    }
    
    # Predictor step
    u<-x + h*omega*t
    q<-qr(jac(u) |> t())  |> qr.Q(complete=TRUE)
    
    disto <- 0
    decel <- 1/max_decel # deceleration factor
    
    for (it in 1:maxiter) {
      y<-H(u)
      
      
      # change in rho proposed by Newton step
      # This is the step we'd use if we kept lambda constant
      #drho <- -solve(Hrho(u[nActions+1],u[1:nActions]),y)
      drho<- - ginv(jac(u)) %*% y
      dist<-(drho^2) |> sqrt() |> sum()
      
      
      if (dist >= max_dist) {
        print(paste("Proposed distance",dist,"exceeds max_dist =",max_dist))
        accept<-FALSE
        break
      }
      
      decel<-max(c(decel,sqrt(dist/max_dist)*max_decel))
      
      if (it>1) {
        contr<- dist/(disto +tol*eta)
        if (contr > max_contr) {
          print(paste("Maximum contraction rate exceeded"))
          accept<-FALSE
          break
        }
        decel <-max(c(decel,sqrt(contr/max_contr)*max_decel))
      }
      if (dist < tol) {
        
        # Success! update and break out of iteration
        break
      }
      disto<-dist
      
      # if we have got to this point, then:
      #  1. The Newton step has not proposed a too-large change
      #  2. The contraction rate has not been exceeded
      #  3. The proposed Newton step is not so small that we have effectively 
      #    found a solution
      # Therefore, we update rho
      
      u<-u+drho
      
      
      
      if (it==maxiter) {
        # We have run out of iterations. Terminate
        print(paste("Maximum number of iterations",maxiter,"reached"))
      }
      
    } ## END OF CORRECTOR STEP
    
    if (!accept) {
      # Step was not accepted; take a smaller step and try again
      h<-h/max_decel
    }
    
    # Standard steplength adaptation
    h<-abs(h/min(c(decel,max_decel)))
    
    # Update with outcome of successful PC step
    if (sum(t*q[,nActions+1])<0) {
      # The orientation of the curve as determined by the QR
      # decomposition has changed.
      #
      # The original Allgower-Georg QR decomposition implementation
      # ensures the sign of the tangent is stable; when it is stable
      # the curve orientation switching is a sign of a bifurcation.
      omega<- -omega
      #print("sign of omega flipped")
    }
    
    
    
    x<-u
    
    # Force the sum-to-one constraint to be a bit more binding
    #x[1:nActions]<-x[1:nActions]-log(sum(exp(x[1:nActions])))
    
    t<-q[,nActions+1]
    
    #print(c(exp(x[1:nActions]),x[nActions+1]))
    
    
    QRE<-rbind(QRE,x |> as.vector())
    #print(c(exp(x[1:nActions]),x[nActions+1]) |> t())
    print(paste("lambda =",x[nActions+1]))
    
    
  }
  QRE
  
  
}

softmax<-function(x) {
  
  x<-x-max(x)
  
  exp(x-log(sum(exp(x))))
  
}

softmaxRHO<-function(x,rho) {
  x<-x-max(x)
  exp(x-rho-log(sum(exp(x))))
}


softmax_vectorized<-function(x) {
  # Takes a matrix input. Each column is a vector of expected payoffs
  # Returns a matrix, each column is a vector of logit responses
  
  x<-x-rep(1,dim(x)[1])%*%(apply(x,MARGIN = 2,FUN=max) |> t())
  
  sumexpx<-apply(exp(x),MARGIN=2,FUN=sum)
  
  exp(x-rep(1,dim(x)[1])%*%t(log(sumexpx)))
  
}

softmaxRHO_vectorized<-function(x,rho) {
  
  x<-x-rep(1,dim(x)[1])%*%(apply(x,MARGIN = 2,FUN=max) |> t())
  
  sumexpx<-apply(exp(x),MARGIN=2,FUN=sum)
  
  exp(x-rho%*%t(rep(1,dim(x)[2]))-rep(1,dim(x)[1])%*%t(log(sumexpx)))
  
}