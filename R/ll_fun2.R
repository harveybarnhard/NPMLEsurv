# Harvey Barnhard
# February 29, 2020
# Last modified on February 29, 2020

# Libraries ====================================================================
library(Rsolnp)   # Constrained optimization
library(parallel) # Parallel processing
library(pbapply)  # Progress bars for parallel processing
library(TMB)

prop.hazard <- function(Xlist,
                        censorvec,
                        thetahat,
                        theta_dom,
                        numiter=8){
  data   <- list(Xlist=Xlist, censorvec=censorvec)

  # Step 1: Initialize coefficient vectors
  alphahat <- rep(0, ncol(Xlist[[1]]))
  pihat    <- 1

  loglik <- c()
  i <- 1
  while(i <= numiter){
    cat("========== Iteration ", i, "==========\n")
    # Set parameters for iteration
    parameters <- list(alpha=alphahat,
                       theta=thetahat,
                       pi=pihat)
    # Step 2: Create objective function and optimize keeping heterogeneity
    # static.
    map <- list(theta=rep(factor(NA), length(thetahat)),
                pi=rep(factor(NA), length(pihat)))
    obj <- TMB::MakeADFun(data,
                          parameters,
                          DLL="NPMLEsurv",
                          map=map,
                          silent=TRUE)
    obj$method  <- "BFGS" # Optimization method
    obj$hessian <- TRUE  # Return Hessian?
    optiter <- do.call("optim", obj)

    # If only one iteration is desired, then return initial values of
    # heterogeneity points
    if(numiter==1){
      break
    }

    # Check to see if the negative log-likelihood has decreased by more than
    # 0.5. If not, end process
    if(i > 1){
      if(abs(optiter$value - loglik[length(loglik)]) < 0.5){
        break
      }
    }
    loglik[i] <- optiter$value
    alphahat  <- optiter$par
    # Print parameter output
    cat(paste0(alphahat, "\n"))

    # Step 3: Evaluate gradient of a new heterogeneity support point over a
    # preset grid of values
    parameters <- list(alpha=alphahat,
                       theta=c(thetahat,0),
                       pi=c(pihat,0))
    map <- list(alpha=rep(factor(NA), length(alphahat)))
    obj <- TMB::MakeADFun(data,
                          parameters,
                          DLL="NPMLEsurv",
                          map=map,
                          silent=TRUE)
    muvec <- sapply(theta_dom[!theta_dom%in%thetahat],
                    function(x) obj$gr(c(thetahat,
                                         x,
                                         c(pihat,0)),
                                       order=1)[2*(length(pihat)+1)])
    if(all(muvec>=0)){
      break
    }
    thetahat <- c(thetahat, theta_dom[!theta_dom%in%thetahat][which.min(muvec)])
    pihat    <- rep(1/(length(thetahat)), length(thetahat))

    # Step 4: Numerically solve the constrained optimization problem for
    # optimal probabilities
    parameters <- list(alpha=alphahat,
                       theta=thetahat,
                       pi=pihat)
    map <- list(alpha=rep(factor(NA), length(alphahat)),
                theta=rep(factor(NA), length(thetahat)))
    obj <- TMB::MakeADFun(data,
                          parameters,
                          DLL="NPMLEsurv",
                          map=map,
                          silent=TRUE)
    eqfun <- function(x) sum(x)
    opt <- solnp(pihat,
                 fun=obj$fn,
                 eqfun = eqfun,
                 eqB =1,
                 LB=rep(0, length(pihat)),
                 UB=rep(1, length(pihat)),
                 control=list(trace=0))
    pihat <- opt$pars
    i <- i+1
  }
  return(list(alpha=alphahat,
              pi=pihat,
              theta=thetahat,
              loglik=loglik,
              fisher=optiter$hessian))
}

# Wrapper ======================================================================
est.prop.hazard <- function(Xlist, censorvec, theta_dom, numiter=8,
                            clust, theta_num){
  # Heterogenity support points to start from
  theta_start <- sample(theta_dom, theta_num, replace=FALSE)

  # Create cluster
  trash <- clusterEvalQ(clust, library("Rsolnp", "TMB"))
  clusterExport(clust,
                c("Xlist", "censorvec", "theta_dom", "theta_start",
                  "prop.hazard", "numiter"),
                envir=environment())
  trash <- clusterEvalQ(clust, dyn.load(TMB::dynlib("NPMLEsurv")))
  # Estimate
  results <- pblapply(theta_start, function (x){
    prop.hazard(
      Xlist=Xlist,
      censorvec,
      x,
      theta_dom,
      numiter
    )
  },
  cl=clust)
  # Find the starting value of theta that resulted in the lowest negative
  # log-likelihood
  loglik  <- sapply(lapply(results, `[[`, 4), min)
  optresults <- results[[which.min(loglik)]]
  # Find the estimates and approximate standard errors using the delta method
  alphahat <- optresults$alpha
  fisher  <- optresults$fisher
  se <- 1/sqrt(diag(fisher))
  # Name coefficients
  names(alphahat) <- paste0("alpha", 1:length(alphahat))
  names(se) <- names(alphahat)
  # Output results
  output <- list(coef=rbind(alphahat, se),
                 ll=min(loglik))
  return(output)
}
