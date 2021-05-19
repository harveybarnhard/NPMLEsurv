# Set paths and libraries ======================================================
basepath <- "C:/Users/Harvey/Dropbox/Personal/UChicago/4_Year/Thesis"
library(data.table) # data manipulation
library(parallel)   # For parallel processing
library(TMB)        # Integrating C++ code
source(paste0(basepath, "/code/clean/inspection/aux_opt/ll_fun2.R"))


# Compile log-likelihood function ==============================================
setwd(paste0(basepath, "/code/clean/inspection/aux_opt"))
compile("loglik.cpp")
dyn.load(dynlib("loglik"))

# Load dataset in list format ==================================================
load(paste0(basepath, "/data/clean/data_lists/t_halfyr_0_0_fail_2_2.Rdata"))
load(paste0(basepath, "/data/clean/data_lists/censorvec.Rdata"))

# Estimate parameters and standard errors ======================================
theta_start <- runif(50, -11, 3) # Grid of heterogeneity points
cl <- makeCluster(12)
results <- est.prop.hazard(Xlist,
                           censorvec,
                           theta_dom=theta_start,
                           clust=cl,
                           theta_num=length(theta_start))
stopCluster(cl)

