#' Give the estimated covariance matrix from the L matrix (see https://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html)
#' @param TMB_obj: TMB object
#' @param name_cov_par: the name of the parameters which correspond to the elements of L
#' @return the covariance matrix
L_to_cov_wrapper_TMB <- function(TMB_obj, name_cov_par='cov_RE', d=NULL){
  .x <- summary(TMB_obj)
  if(is.null(d)){
    .d <- nrow(python_like_select_rownames(.x, 'beta'))/2
  }else{
    .d <- d
  }
  L_to_cov(python_like_select_rownames(.x, name_cov_par)[,1], d = .d)
}

#' @param object: TMB-input list object
#' @return  TMB-input list object with columns sorted by increasing colsum
sort_columns_TMB = function(object){
  order_cats <- order(colSums(object$Y), decreasing = F) ## increasing so that the last category is the highest
  object$Y = object$Y[,order_cats]
  return(object)
}

#' @param object: TMB-input list object
#' @return  TMB-input list object with columns sorted such that SBS1 is the category in the last column
sort_columns_TMB_SBS1 = function(object){
  if("SBS1" %in% colnames(object$Y)){
    order_cats <- c(which('SBS1' != colnames(object$Y)), which('SBS1' == colnames(object$Y)))
    object$Y = object$Y[,order_cats]
  }else{
    warning('There are no exposures for SBS1 in this sample. Keeping the same order')
  }
  return(object)
}

#' Function to run the model using TMB
#' @param model: name of model
#' @param object: object of type exposures_inputTMB (including X, Y, Z)
#' @param smart_init_vals: boolean, whether a fixed-effects multinomial regression should be run first to get initial estimates
#' @param use_nlminb: boolean, whether <nlminb> should be used for estimation. Alternatively, <optim> is used
#' @param initial_params: (optional) list of initial parameters for estimation
#' @useDynLib diagRE_dirichletmultinomial_single_lambda
#' @useDynLib diagRE_ME_dirichletmultinomial
#' @useDynLib diagRE_ME_multinomial
#' @useDynLib FE_dirichletmultinomial_single_lambda
#' @useDynLib FE_dirichletmultinomial
#' @useDynLib fullRE_dirichletmultinomial_single_lambda
#' @useDynLib fullRE_ME_dirichletmultinomial_onefixedlambda
#' @useDynLib fullRE_ME_dirichletmultinomial
#' @useDynLib fullRE_ME_multinomial
#' @useDynLib singleRE_dirichlet_multinomial
#' @useDynLib functions.hpp
#' @importFrom Rcpp evalCpp
wrapper_run_TMB = function(model, object=NULL, smart_init_vals=T, use_nlminb=F, initial_params=NULL){
  ## sort_columns=F, 
  
  ## if the object of data and covariates is an argument
  data = object
  
  # if(sort_columns){
  #   order_cats = order(colSums(data$Y), decreasing = T)
  #   data$Y = data$Y[,order_cats]
  # }
  
  data$Y = matrix(data$Y, nrow=nrow(data$Y))
  data$x = (matrix(data$x, ncol=2))
  
  d <- ncol(data$Y) ## number of signatures
  n <- ncol(data$z) ## number of INDIVIDUALS, not samples
  
  if(smart_init_vals){
    require(nnet)
    .x_multinom = multinom(data$Y ~ data$x[,2])
    beta_init = t(coef(.x_multinom))
  }else{
    beta_init = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
                        nrow = 2, byrow=TRUE))
  }
  
  parameters <- list(
    beta = beta_init,
    u_large = matrix(rep(1, (d-1)*n), nrow=n),
    logs_sd_RE=rep(1, d-1),
    cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2)
  )
  
  if(model == "fullRE_M"){
    data$num_individuals = n
    dll_name <- "fullRE_ME_multinomial"
    rdm_vec <- "u_large"
  }else if(model == "diagRE_M"){
    data$num_individuals = n
    dll_name <- "diagRE_ME_multinomial"
    rdm_vec <- "u_large"
  }else if(model == "FE_DM"){
    # data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    data$lambda_accessory_mat = cbind(data$x[,2], 1-data$x[,2])
    
    parameters <- c(parameters, list(log_lambda = matrix(c(2,2))))
    parameters$u_large = NULL
    parameters$logs_sd_RE = NULL
    parameters$cov_par_RE = NULL
    dll_name <- "FE_dirichletmultinomial"
    rdm_vec <- NULL
  }else if(model == "fullRE_DM"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    parameters <- c(parameters,list(log_lambda = matrix(c(2,2))))
    dll_name <- "fullRE_ME_dirichletmultinomial"
    rdm_vec <- "u_large"
  }else if(model == "fullRE_DMonefixedlambda"){
    ## fixing one of the overdispersion parameters
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    parameters <- c(parameters,list(log_lambda = matrix(c(2))))
    dll_name <- "fullRE_ME_dirichletmultinomial_onefixedlambda"
    rdm_vec <- "u_large"
  }else if(model == "diagRE_DM"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    parameters$log_lambda = matrix(c(2,2))
    parameters$cov_par_RE = NULL
    dll_name <- "diagRE_ME_dirichletmultinomial"
    rdm_vec <- "u_large"
  }else if(model  == "fullREDMsinglelambda"){
    data$num_individuals = n
    parameters$log_lambda = 1.1
    rdm_vec <- "u_large"
    dll_name <- "fullRE_dirichletmultinomial_single_lambda"
  }else if(model  == "singleRE_DM"){
    ## single value for RE, shared for all log-ratios
    cat('single value for RE, shared for all log-ratios\n')
    data$num_individuals = n
    print(parameters)
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    parameters$u_random_effects = as(parameters$u_large[,1], 'matrix')
    parameters$u_large = NULL
    parameters$logs_sd_RE = NULL
    parameters$cov_par_RE = NULL
    parameters$logSigma_RE <- 1.1
    parameters <- c(parameters,list(log_lambda = matrix(c(2,2))))
    rdm_vec <- "u_random_effects"
    dll_name <- "singleRE_dirichlet_multinomial"
  }else if(model == "fullREhalfDM"){
    stop("fullREhalfDM used to be done with  fullRE_dirichletmultinomial_single_lambda")
    data$num_individuals = n
    
    parameters <- list(parameters, log_lambda = 1.1)
    rdm_vec <- "u_large"
    dll_name <- "CHANGETHIS"
  }else if(model == "diagREDMsinglelambda"){
    data$num_individuals = n
    data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
    
    parameters <- c(parameters, log_lambda = 1.1)
    parameters$cov_par_RE = NULL
    
    # Error below indicates   <parameters> list not being created correcty
    # Error in MakeADFun(data = Data, parameters = Parameters, random = Random)
    # :
    #   Only numeric matrices, vectors, arrays or factors can be interfaced
    # 
    # 
    # parameters <- list(
    #   beta = beta_init,
    #   u_large = matrix(rep(1, (d-1)*n), nrow=n),
    #   logs_sd_RE=rep(1, d-1),
    #   log_lambda = 2
    # )
    dll_name <- "diagRE_dirichletmultinomial_single_lambda"
    rdm_vec <- "u_large"
    
    # obj <- MakeADFun(data, parameters, DLL="", random = "")
  }else if(model == "FEDMsinglelambda"){
    data$num_individuals = NULL
    parameters$u_large = NULL
    parameters$logs_sd_RE = NULL
    parameters$cov_par_RE = NULL
    parameters <- c(parameters, list(log_lambda = 1.1))
    
    dll_name <- "FE_dirichletmultinomial_single_lambda"
    rdm_vec <- NULL
    
  }else if(model == "fullRE_multinomial_REv2"){
    data$num_individuals = n
    parameters$u_large = NULL
    parameters <- list(parameters,
                       u_large1 = matrix(runif(n)),
                       u_large2 = matrix(runif(n)),
                       u_large3 = matrix(runif(n)),
                       u_large4 = matrix(runif(n)))
    
    dll_name <- "fullRE_ME_multinomial_REv2"
    rdm_vec <- c("u_large1", "u_large2","u_large3", "u_large4")
  }else{
    stop('Specify correct <model>\n')
  }
  
  if(!is.null(initial_params)){
    ## initial parameters are passed as arguments
    parameters <- initial_params
  }
  
  require(TMB)
  # location_package <- find.package('CompSign')
  # TMB::compile(paste0(location_package, "models/", dll_name),  "-std=gnu++17")
  # dyn.load(dynlib(gsub(".cpp", "", paste0(location_package, "models/", dll_name))))
  
  if(is.null(rdm_vec)){
    ## fixed effects model
    obj <- MakeADFun(data, parameters, DLL=dll_name)
  }else{
    ## random effects model
    obj <- MakeADFun(data, parameters, DLL=dll_name, random = rdm_vec)
  }
  
  if(use_nlminb){
    opt = nlminb(start = obj$par, obj = obj$fn, gr = obj$gr, iter.max=iter.max)
  }else{
    obj$hessian <- TRUE
    opt <- do.call("optim", obj)
    opt
    opt$hessian ## <-- FD hessian from optim
  }
  return_report <- sdreport(obj)
  
  return(return_report)
}


wald_generalised = function(v, sigma){
  # warning('20201218: sigma**(1/2) has now been replaced by (as we had before sometime in November) sigma')
  chisqrt_stat = t(v) %*% solve(sigma) %*% v
  pchisq(q = chisqrt_stat, df = length(v), lower.tail = FALSE)
}

#' Generalised wald test for the slope beta values
#' @param i: TMB output object
#' @param verbatim: verbatim boolean
#' @param fail_non_converged: should a non-converged run lead to an NA value?
wald_TMB_wrapper = function(i, verbatim=TRUE, fail_non_converged=T){
  if(typeof(i) == "character"){
    return(NA)
  }else{
    idx_beta = select_slope_2(which(names(i$par.fixed) == "beta"), verbatim=verbatim)
    if(!i$pdHess & fail_non_converged){
      ## didn't converge
      NA
    }else{
      if(length(idx_beta) == 1){
        if(is.na(idx_beta)){
          NA
        }else{
          if(verbatim) cat('Check data - slope appears to be of length one (binomial)\n')
          wald_generalised(v = i$par.fixed[idx_beta], sigma = i$cov.fixed[idx_beta,idx_beta])
        }
      }else{
        wald_generalised(v = i$par.fixed[idx_beta], sigma = i$cov.fixed[idx_beta,idx_beta])
      }
    }
  }
}

load_models <- function(i){
  require(TMB)
  TMB::compile("../R/mm_multinomial/fullRE_ME_dirichletmultinomial.cpp",  "-std=gnu++17")
  dyn.load(dynlib("../R/mm_multinomial/fullRE_ME_dirichletmultinomial"))
}