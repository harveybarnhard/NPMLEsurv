library(TMB)

# Monte-Carlo testing ==========================================================
n     <- 10000
lambda <- 0.05 # Scale parameter
v      <-  2  # Shape parameter
beta1   <- -0.2
beta2  <- 0.3

for(i in 1:1000){
  # Simulate observed and unobserved heterogenetiy parameters
  x1  <- c(rep(0, n/4), rep(1, n/4), rep(0, n/4), rep(1, n/4))
  x2  <- c(rep(0, n/4), rep(0, n/4), rep(1, n/4), rep(1, n/4))
  thetavec <- sample(c(-0.2,0.2), n, replace=TRUE)
  #t  <- ceiling((-log(runif(n))/(lambda*exp(beta1*x1 + beta2*x2 + thetavec)))^(1/v))
  t0 <- sample(2:20, n, replace=TRUE)
  u  <- runif(n, 0, 1)
  t <- c()
  for(j in 1:n){
    if(-log(u[j]) < lambda*exp(beta1*x1[j] + thetavec[j])*t0[j]^v){
      t[j]  <- ceiling((-log(u[j])/(lambda*exp(beta1*x1[j] + thetavec[j])))^(1/v))
    }else{
      t[j]  <- ceiling(((-log(u[j]) -
                            lambda*exp(beta1*x1[j] + thetavec[j])*t0[j]^v +
                            lambda*exp(beta2)*exp(beta1*x1[j] + thetavec[j])*t0[j]^v
                          )/(lambda*exp(beta2)*exp(beta1*x1[j] + thetavec[j])))^(1/v))
    }
  }
  Xlist <- list()
  censorvec <- rep(NA, n)
  censor_time <- 10
  for(j in 1:n){
      # Xlist[[j]] <- cbind(rep(x1[j], t[j]),
      #                     as.integer((1:t[j]) >= floor(t0[j])),
      #                     rep(1, t[j]),
      #                     cbind(diag(1, nrow=t[j], ncol=max(t))))
      Xlist[[j]] <- cbind(rep(x1[j], t[j]),
                          as.integer((1:t[j]) >= floor(t0[j])),
                          rep(1, t[j]),
                          0:(t[j]-1))
      censorvec[j] <- t[j]>=censor_time
  }
  theta_dom <- runif(50, min=-5, max=3)
  thetahat_init  <- runif(100, min=-5, max=5)
# Set parameter bounds =========================================================
  # Set bounds for all parameter
  L = c(alpha=3,theta=0,pi=1)
  U = c(alpha=5,theta=0, pi=1)

  # List parameters that should be fixed in nlminb() or optim()
  map = list(theta=factor(NA), pi=factor(NA))
  #map=list()  # All parameters active

  # Remove inactive parameters from bounds
  (L <- L[-match(names(map), names(L))])
  (U <- U[-match(names(map), names(U))])
# Fit the model ================================================================
  data <- list(Xlist=Xlist, censorvec=censorvec)
  parameters <- list(alpha=c(0,0,0,0), theta=0, pi=1)
  obj <- TMB::MakeADFun(data, parameters, DLL="NPMLEsurv", map=map)
  obj$method  <- "BFGS"
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  parframe <- rbind(parframe, opt$par)
}
