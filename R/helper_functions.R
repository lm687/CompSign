#' Compute the softmax-transformed version of a vector or a matrix. Reminder: if x is an ALR-transformed compositional vector, you want to append 0 as the final element of the vector (if x is a vector), or a final column of zeros (if x is a matrix)
#' @param x: a vector or a matrix
#' @return the softmax-transformed x
softmax = function(x){
  if(is.null(dim(x))){
    ## vector
    sum_x = sum(exp(x))
    exp(x)/sum_x
  }else{
    ## matrix
    sum_x = rowSums(exp(x))
    sweep(exp(x), 1, sum_x, '/')
  }
}

split_rows <- function(df, f){
  lapply(unique(f), function(fact) df[which(f == fact),])
}

reflect_matrix = function(m){
  m[nrow(m):1,]
}

#' Give a Z matrix for a scenario in which each patient contributes to two samples, and the samples are sorted by group, then by patient.
#' @param n_times_2: the number of samples times two
#' @return the Z matrix
give_z_matrix = function(n_times_2){
  a = matrix(0, nrow = n_times_2/2, ncol = n_times_2/2)
  diag(a) = 1
  rbind(a, a)
}

#' Give a X (covariate) matrix for a scenario in which each patient contributes to two samples, and the samples are sorted by group, then by patient.
#' @param n_times_2: the number of samples times two
#' @return the X matrix
give_x_matrix = function(n_times_2){
  cbind(rep(1, n_times_2), rep(c(0,1), each=n_times_2/2))
}

rsq = function (x, y) cor(x, y) ^ 2


python_like_select = function(vector, grep_substring){
  vector[grepl(pattern = grep_substring, x = vector)]
}

python_like_select_name = function(vector, grep_substring){
  vector[grepl(pattern = grep_substring, x = names(vector))]
}

python_like_select_colnames = function(matrix, grep_substring){
  matrix[,grepl(pattern = grep_substring, x = colnames(matrix))]
}

python_like_select_rownames = function(matrix, grep_substring){
  matrix[grepl(pattern = grep_substring, x = rownames(matrix)),]
}

clean_name = function(x){
  gsub(".RDS", "", paste0(strsplit(x, "_")[[1]][2:3], collapse = ""))
}

clean_name_fullRE = function(x){
  gsub(".RDS", "", paste0(strsplit(x, "_")[[1]][3:4], collapse = ""))
}

clean_name_fullRE_2 = function(x){
  gsub(".RDS", "", paste0(strsplit(x, "_")[[1]][4:5], collapse = ""))
}

re_vector_to_matrix = function(vec_RE, dmin1){
  ## Random effects vector to matrix
  matrix(vec_RE, ncol=dmin1)
}

#' @param object: TMB result object
#' @return  summary indicating if the model has converged
give_summary_per_sample = function(TMB_object, verbatim=T){
  if(is.null(TMB_object)){
    "Object doesn't exist"
  }else{
    if(typeof(TMB_object) %in% c("character", "logical")){
      return('Timeout or some error')
    }else{
      if(TMB_object$pdHess){
        return('Good')
      }else{
        return('Non-PD')
      }
    }
  }
}

#' @param object: list of TMB results objects
#' @return  list summary indicating if the models have converged
give_summary_of_runs2 = function(vector_TMB_objects, long_return, verbatim=T){
  timeout_bool = sapply(vector_TMB_objects, typeof) == "character"
  hessian_positivedefinite_bool = sapply(vector_TMB_objects[!timeout_bool], function(i){
    if(length(i) == 1){FALSE}else{i$pdHess}})
  summary_runs = c(sum(timeout_bool), sum(!hessian_positivedefinite_bool), sum(hessian_positivedefinite_bool))
  names(summary_runs) = c('(failed) timeout', '(failed) non-positive semi-definite hessian', '(successful) positive semi-definite hessian')
  if(long_return){
    list(Timeout=names(vector_TMB_objects)[timeout_bool],
         hessian_nonpositivedefinite_bool = names(vector_TMB_objects)[!timeout_bool][!hessian_positivedefinite_bool],
         hessian_positivedefinite_bool = names(vector_TMB_objects)[!timeout_bool][hessian_positivedefinite_bool])
  }else{
    return(summary_runs)
  }
}


#' @param object: TMB result object
#' @return  p-value for wald test for overdispersion (h1= lambda \neq 0)
wald_TMB_wrapper_overdisp = function(i, verbatim=TRUE){
  ## wald test for the overdispersion parameter
  if(typeof(i) == "character"){
    return(NA)
  }else{
    idx_beta = which(names(i$par.fixed) == "log_lambda")
    if(!i$pdHess){
      ## didn't converge
      NA
    }else{
      scaled_coefs <- scale(i$par.fixed[idx_beta], center = T, scale = F)
      wald_generalised(v = as.vector(scaled_coefs),
                       sigma = i$cov.fixed[idx_beta,idx_beta])
    }
  }
}

#' @param object: TMB result object
#' @return  p-value for t-test for overdispersion (h1= lambda \neq 0)
ttest_TMB_wrapper_overdisp = function(i, verbatim=TRUE){
  ## wald test for the overdispersion parameter
  if(typeof(i) == "character"){
    return(NA)
  }else{
    idx_beta = which(names(i$par.fixed) == "log_lambda")
    if(!i$pdHess){
      ## didn't converge
      NA
    }else{
      .summary <- summary(i)
      loglambdas <- python_like_select_rownames(.summary, 'log_lambda')
      mean_coefs <- loglambdas[,1]
      SE_coefs <- loglambdas[,2]
      
      dmin1 <- nrow(python_like_select_rownames(.summary, 'beta'))/2
      num_patients <- nrow(python_like_select_rownames(.summary, 'u_large'))/dmin1
      tstatistic <- (mean_coefs[1]-mean_coefs[2])/sum(SE_coefs)
      df <- num_patients*2 - 2
      #' For significance testing, the degrees of freedom for this test
      #' is 2n − 2 where n is the number of participants in each group. 
      2*pt(-abs(tstatistic), df=df)
    }
  }
}

#' @param object: TMB result object
#' @return  subset of beta estimates for the slope
select_slope_2 = function(i, verbatim=TRUE){
  if(is.null(dim(i))){
    # if(verbatim) warning('As per 27 August 2020 it seems clear that this version, and not <select_slope>, is correct')
    i[c(F,T)]
  }else{
    i[,c(F,T)]
  }
}

#' @param object: TMB result object
#' @return  subset of beta estimates for the intercept
select_intercept = function(i, verbatim=TRUE){
  if(is.null(dim(i))){
    i[c(T,F)]
  }else{
    i[,c(T,F)]
  }
}


vector_to_ct_list = function(vec){
  ##' given a vector which contains two or more types of features for the same cancer types,
  ##' convert it to a matrix with the pairing per cancer type
  which_sigs = (grep('signatures', names(vec)))
  which_nuc = (grep('nucleotidesubstitution1', names(vec)))
  if(sum(length(which_sigs)+length(which_nuc)) < length(vec)){stop('There is a third feature category')}
  ct_naked1 = gsub('signatures', '', names(vec)[which_sigs])
  ct_naked2 = gsub('nucleotidesubstitution1', '', names(vec)[which_nuc])
  ct_naked1[match(ct_naked1, ct_naked2)]
  ret = cbind(vec[which_sigs][match(ct_naked1, ct_naked2)], vec[which_nuc])
  rownames(ret) = gsub('signatures', '', rownames(ret))
  colnames(ret) = c('signatures', 'nucleotidesubsitution1')
  return(ret)
}

give_UNSTRUCTURED_CORR_t_matrix = function(vec, dim_mat){
  # #https://kaskr.github.io/adcomp/classUNSTRUCTURED__CORR__t.html
  m = matrix(1, nrow = dim_mat, ncol = dim_mat)
  ## fill in the order that TMB's UNSTRUCTURED_CORR_t saves the covariances
  m[unlist(sapply(2:nrow(m), function(rw) seq(from = rw,length.out = (rw-1), by = nrow(m) )))] = vec
  m[unlist(sapply(2:nrow(m), function(cl) seq(from = (((cl-1)*nrow(m))+1),length.out = (cl-1), by = 1 )))] = vec
  return(m)
}

get_count_object = function(ct, feature, pre_path=NULL){
  if(is.null(pre_path)){
    pre_path = "../../data/roo/"
  }
  readRDS(paste0(pre_path, ct, "_", feature, "_ROO.RDS"))
}

get_count_object_file = function(fle){
  readRDS(fle)
}

#' @param object: TMB result object
#' @return standard errors of beta estimates
give_stderr = function(i, only_slopes=T, only_betas=T){
  if(!only_betas){
    stop('Not yet implemented')
  }else{
    if(length(i) == 1){ ## they are just NAs
      if(only_slopes){
        .x = select_slope_2(python_like_select_name(i$par.fixed, "beta"), verbatim = FALSE) ## repeat the NAs
      }else{
        .x = python_like_select_name(i$par.fixed, "beta") ## repeat the NAs
      }
    }else{
      if(only_slopes){
        .x = select_slope_2(python_like_select_name(TMB::summary.sdreport(i)[,2], "beta"),v=F)
      }else{
        .x = python_like_select_name(TMB::summary.sdreport(i)[,2], "beta")
      }
    }
    .x
  }
}

createbarplot_object = function(fle, slotname='count_matrices_all', pre_path_funs=NULL){
  .x = readRDS(fle)
  if(is.null(pre_path_funs)) pre_path_funs <- "../../../CDA_in_Cancer/code/"
  lapply(.x, createbarplot_ROOSigs, slot=slotname, pre_path = pre_path_funs)
}


simulate_from_DM_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, integer_overdispersion_param){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  overdispersion_lambda = integer_overdispersion_param*exp(python_like_select_name(tmb_fit_object$par.fixed, "log_lambda")[x_matrix[,2]+1])
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* softmax(c(logRmat[l,], 0)))))
  }else{
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* sim_thetas[l,])))
  }
  return(sim_thetas)
}

give_dummy_colnames <- function(i){colnames(i) <- paste0('c', 1:ncol(i)); i}
give_dummy_rownames <- function(i){rownames(i) <- paste0('r', 1:nrow(i)); i}

resort_columns <- function(i, order){
  i$Y = i$Y[,order]
  i
}

simulate_from_M_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = softmax(cbind(logRmat, 0))
  }else{
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
  }
  return(sim_thetas)
}

simulate_from_correlated_binom = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, return_logratios=F){
  d = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, d)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    if(return_logratios){
      return(logRmat)
    }else{
      ## return first probability
      return(apply(logRmat, 2, function(i) exp(i)/(1+exp(i))))
    }
  }else{
    stop('Not implemented yet')
    # sim_thetas = sapply(1:d,
    #                     function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
    #                              give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    return(sim_thetas)
  }
}

#' @param object: vector of estimantes and vector of standard errors
#' @return confidence interval for each estimate 
give_confidence_interval = function(vec_est, vec_stderr){
  sapply(1:length(vec_est), function(i) c(vec_est[i]-1.96*vec_stderr[i],vec_est[i]+1.96*vec_stderr[i]) )
}

#' @param object: vector of estimantes and vector of standard errors, and another vector for testing
#' @return boolean testing whether each element of <vec_true> is included in the confidence interval for each estimate
give_params_in_CI = function(vec_est, vec_stderr, vec_true){
  ci = give_confidence_interval(vec_est, vec_stderr)
  return(sapply(1:length(vec_est), function(i) (ci[1,i] < vec_true[i]) & (ci[2,i] > vec_true[i]) ))
}

wrapper_run_ttest_ilr = function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  props = sapply(x, normalise_rw, simplify = FALSE)
  return(Compositional::hotel2T2(x1 = compositions::ilr(props[[1]]), x2 = compositions::ilr(props[[2]]))$pvalue)
}
wrapper_run_ttest_props = function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  props = sapply(x, normalise_rw, simplify = FALSE)
  return(Compositional::hotel2T2(x1 =props[[1]][,-1], x2 = props[[2]][,-1])$pvalue)
}


summarise_DA_detection = function(true, predicted, verbose=T){
  if(verbose) warning('11 August 2021: there was a problem in how FP, TP, etc. were calculated')
  require(ROCR)
  ## remove NAs
  which_na = which(is.na(predicted))
  if(length(which_na) > 0){ ## some NA
    true = true[-which_na]
    predicted = predicted[-which_na]
  }
  
  FPs = sum(!true & predicted)
  FPR = FPs/sum(!true)
  TPs = sum(true & predicted)
  TPR = TPs/sum(true)
  TNs = sum(!true & !predicted)
  TNR = TNs/sum(!true)
  FNs = sum(true & !predicted)
  FNR = FNs/sum(!true)
  Accuracy=(TPs + TNs)/length(true)
  if(sum(true) == 0 | sum(!true) == 0){
    WeightedAccuracy = NA
  }else{
    WeightedAccuracy=mean(c(TPs/sum(true), TNs/sum(!true)))
  }
  total_pos = sum(true | predicted)
  Power = TPs/total_pos
  Sensitivity = TPs / (TPs + FNs)
  Specificity = TNs / (TNs + FPs)
  Recall=Sensitivity
  Precision=TPs/(TPs + FPs)
  pred <- (try(ROCR::prediction(as.numeric(true), as.numeric(predicted))))
  if(typeof(pred) == 'S4'){
    AUC = as.numeric(try(ROCR::performance(pred, "auc")@y.values[[1]]))
  }else{
    AUC = NA
  }
  return(c(FPR=FPR, TPR=TPR, Power=Power, AUC=AUC, Specificity=Specificity,
           Sensitivity=Sensitivity, Recall=Recall, Precision=Precision,
           Accuracy=Accuracy, WeightedAccuracy=WeightedAccuracy))
}

fill_covariance_matrix = function(arg_d, arg_entries_var, arg_entries_cov, verbose=T){
  if(verbose) warning('This function had been incorrect until now (30 july 2021)')
  .sigma <- give_UNSTRUCTURED_CORR_t_matrix(vec = arg_entries_cov, dim_mat = arg_d)
  # .sigma = matrix(NA, arg_d, arg_d)
  diag(.sigma) = arg_entries_var
  # .sigma[unlist(sapply(1:(arg_d-1), function(i) (i-1)*arg_d + (i+1):arg_d ))] = arg_entries_cov
  # .sigma[unlist(sapply(1:(arg_d-1), function(i) (i) + ((i):(arg_d-1))*arg_d))] = arg_entries_cov
  return(.sigma)
}

extract_sigs_TMB_obj <- function(dataset_obj_trinucleotide, subset_signatures, signature_version='v3', signature_definition=NULL){
  require(mutSigExtractor)
  
  if(!is.null(signature_definition)){
    cat('Using specified signature definitions\n')
    sigdefs <- signature_definition
  }else{
    if(signature_version == 'v3'){
      sigdefs <- SBS_SIGNATURE_PROFILES_V3
    }else if(signature_version == 'v2'){
      sigdefs <- SBS_SIGNATURE_PROFILES_V2
      colnames(sigdefs) <- gsub("SBS0", "SBS", colnames(sigdefs))
      print(colnames(sigdefs))
    }else{
      stop('Incorrect version of signatures')
    }
  }
  
  if(length(subset_signatures) == 1){
    sigdefs_subset <- select(sigdefs, subset_signatures)
  }else{
    sigdefs_subset <- sigdefs[,match(subset_signatures, colnames(sigdefs))]
  }
  dataset_obj_trinucleotide
  
  if(length(subset_signatures) > 1){
    sigdefs_subset = sigdefs_subset[match(colnames(dataset_obj_trinucleotide$Y), rownames(sigdefs_subset)),]
  }else{
    sigdefs_subset = as(sigdefs_subset[match(colnames(dataset_obj_trinucleotide$Y), rownames(sigdefs_subset)),], 'matrix')
    colnames(sigdefs_subset) <- subset_signatures
    rownames(sigdefs_subset) <- colnames(dataset_obj_trinucleotide$Y)
  }
  
  sigs <- fitToSignatures(
    mut.context.counts=(dataset_obj_trinucleotide$Y), 
    signature.profiles=sigdefs_subset
  )
  
  dataset_obj_new <- dataset_obj_trinucleotide
  dataset_obj_new$Y = round(sigs)
  return(dataset_obj_new)
}

compare_signaturefit_to_data <- function(tmb_obj_exposures, tmb_obj_trinucleotide, signature_defs, only_cosim=F, plot=F){
  require(lsa)
  signature_defs
  signature_defs <- signature_defs[match( colnames(tmb_obj_trinucleotide$Y), rownames(signature_defs)),]
  
  reconstructed_trinuc <- as(tmb_obj_exposures$Y, 'matrix') %*% t(as(signature_defs[,colnames(tmb_obj_exposures$Y)], 'matrix'))
  
  if(plot){
    plot(unlist(tmb_obj_trinucleotide$Y), unlist(reconstructed_trinuc))
  }
  
  if(only_cosim){
    cossim <- ( lsa::cosine(x = normalise_rw(colSums(tmb_obj_trinucleotide$Y)), y = normalise_rw(colSums(reconstructed_trinuc))))
    return(cossim)
  }else{
    
    rss <- sum( (tmb_obj_trinucleotide$Y - reconstructed_trinuc)**2)
    rss_max <- max(rowSums( (tmb_obj_trinucleotide$Y - reconstructed_trinuc)**2))
    rss_norm <- sum( (normalise_rw(tmb_obj_trinucleotide$Y) - normalise_rw(reconstructed_trinuc))**2)
    rss_norm_max <- max(rowSums( (normalise_rw(tmb_obj_trinucleotide$Y) - normalise_rw(reconstructed_trinuc))**2))
    cossim <- ( lsa::cosine(x = normalise_rw(colSums(tmb_obj_trinucleotide$Y)), y = normalise_rw(colSums(reconstructed_trinuc))))
    
    return(list(rss=rss, cossim=cossim, rss_max=rss_max, rss_norm=rss_norm, rss_norm_max=rss_norm_max))
  }
}

give_plot_fits_wrapper <- function(list_tmb_obj_exposures,  tmb_obj_trinucleotide, signature_defs){
  .fits_sigs <- sapply(list_tmb_obj_exposures, function(it){
    compare_signaturefit_to_data(it, tmb_obj_trinucleotide, signature_defs)
  })
  colnames(.fits_sigs) <- lapply(list_tmb_obj_exposures, function(i) paste(colnames(i$Y), collapse = '_'))
  .a <- t(.fits_sigs)
  
  .a <- data.frame(RSS=unlist(.a[,'rss']), CosSim=unlist(.a[,'cossim']), rss_max=unlist(.a[,'rss_max']),
                   rss_norm=unlist(.a[,'rss_norm']), rss_norm_max=unlist(.a[,'rss_norm_max']))
  nm <- gsub("SBS", "", rownames(.a))
  .a <- (melt(.a))
  .a$sig <- nm
  .a$sig <- factor(.a$sig, levels=nm)
  ggplot(.a, aes(x=sig, y=value, col=sig))+facet_wrap(.~variable, scales="free_y")+geom_point()+
    geom_line(aes(group=1), col='black')+theme_bw()+
    theme(legend.position = "bottom")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    #guides(col=guide_legend(ncol=1,byrow=TRUE))+
    labs(col='Signature subset')
}

give_ggplot_sig_cor_TMB <- function(tmb_obj_1, tmb_obj_2, ct_it='', no_guides=F){
  mat1 <- tmb_obj_1$Y
  mat2 <- tmb_obj_2$Y
  mat2 <- mat2[,remove_na(match(colnames(mat1), colnames(mat2)))]
  mat1 <- mat1[,remove_na(match(colnames(mat2), colnames(mat1)))]
  
  if(! not_same_sigs){
    if(nrow(mat1) != nrow(mat2)){
      return(ggplot()+ggtitle('Not the same number of samples\n'))
    }
    title=paste0('All sigs\n', ct_it, '\n Not the same number of sigs')
  }else{
    title <- paste0('All sigs\n', ct_it)
  }
  
  df_MSE_comparison <- data.frame(MSE=as.vector(mat1), QP=as.vector(mat2),
                                  sig=rep(colnames(mat1), each=nrow(mat1)))
  plt <- ggplot(df_MSE_comparison, aes(x=MSE, y=QP, col=sig))+
    geom_abline(intercept = 0, slope = 1, lty='dashed')+
    geom_point()+theme_bw()+ggtitle(title)+labs(x='Exposures method 1', y='Exposures method 2', col='Signatures')
  if(no_guides) plt+ guides(col=F)
  print(plt)
}

#' @param sig_obj: signatures object
#' @return signatures object with certain categories (signatures) removed
give_subset_sigs = function(sig_obj, sigs_to_remove){
  
  if(typedata %in% c("nucleotidesubstitution1", "nucleotidesubstitution3", "simulation")){
    slot_name = "count_matrices_all"
  }else if(grepl("signatures", typedata)){
    if(is.null(attr(sig_obj,"count_matrices_active")[[1]]) | (length(attr(sig_obj,"count_matrices_active")[[1]]) == 0)){
      ## no active signatures
      slot_name = "count_matrices_all"
    }else{
      slot_name = "count_matrices_active"
    }
  }
  sig_obj_slot = attr(sig_obj,slot_name)
  sig_obj_slot = lapply(sig_obj_slot, function(i) i[, !(colnames(i) %in% sigs_to_remove)])
  # sig_obj$Y = sig_obj$Y[, !(colnames(sig_obj$Y) %in% sigs_to_remove)]
  attr(sig_obj,slot_name) = sig_obj_slot
  return(sig_obj)
}

#' @param sig_obj: signatures object
#' @return signatures object with certain categories (signatures) amalgamated
give_amalgamated_exposures = function(sig_obj, list_groupings){
  
  if(typedata %in% c("nucleotidesubstitution1", "nucleotidesubstitution3", "simulation")){
    slot_name = "count_matrices_all"
  }else if(grepl("signatures", typedata)){
    if(is.null(attr(sig_obj,"count_matrices_active")[[1]]) | (length(attr(sig_obj,"count_matrices_active")[[1]]) == 0)){
      ## no active signatures
      slot_name = "count_matrices_all"
    }else{
      slot_name = "count_matrices_active"
    }
  }
  sig_obj_slot = attr(sig_obj,slot_name)
  sig_obj_slot = lapply(sig_obj_slot, function(i){
    new_mat = sapply(list_groupings, function(j){
      grouped_exp <- i[,colnames(i) %in% j]
      if(!is.null(ncol(grouped_exp))){
        rowSums(grouped_exp)
      }else{
        grouped_exp
      }
    })
    colnames(new_mat) = sapply(list_groupings, function(i){if(length(i)>1){paste0(i[1], '+')}else{i}})
    return(new_mat)
  })
  attr(sig_obj,slot_name) = sig_obj_slot
  return(sig_obj)
}

#' @param sig_obj: TMB input list object
#' @return TMB input list object with certain categories (signatures) amalgamated
give_amalgamated_exposures_TMBobj = function(sig_obj, list_groupings){
  sig_obj$Y = sapply(list_groupings, function(j){
    grouped_exp <- sig_obj$Y[,colnames(sig_obj$Y) %in% j]
    if(!is.null(ncol(grouped_exp))){
      rowSums(grouped_exp)
    }else{
      grouped_exp
    }
  })
  colnames(sig_obj$Y) = sapply(list_groupings, function(i){if(length(i)>1){paste0(i[1], '+')}else{i}})
  return(sig_obj)
}


give_subset_sigs_TMBobj = function(sig_obj, sigs_to_remove, remove_zero_rows=T){
  sig_obj$Y = sig_obj$Y[,!(colnames(sig_obj$Y) %in% sigs_to_remove)]
  if(remove_zero_rows){
    keep_obs = rowSums(sig_obj$Y) > 0
  }else{
    keep_obs <- 1:nrow(sig_obj$Y)
  }
  sig_obj$Y = sig_obj$Y[keep_obs,]
  sig_obj$d = ncol(sig_obj$Y) ## 20220222
  sig_obj$x = sig_obj$x[keep_obs,]
  sig_obj$z = sig_obj$z[keep_obs,]
  return(sig_obj)
}

give_subset_samples_TMBobj = function(sig_obj, samples_to_remove){
  if(is.null(rownames(sig_obj$Y))){
    keep_samps <- (1:nrow(sig_obj$Y))[-samples_to_remove]
  }else{
    keep_samps <- !(rownames(sig_obj$Y) %in% samples_to_remove)
  }
  sig_obj$Y = sig_obj$Y[keep_samps,]
  sig_obj$x = sig_obj$x[keep_samps,]
  sig_obj$z = sig_obj$z[keep_samps,]
  sig_obj$z <- sig_obj$z[,colSums(sig_obj$z)>0]
  sig_obj$num_individuals <- ncol(sig_obj$z)
  sig_obj$n <- nrow(sig_obj$z) ## 20220222
  return(sig_obj)
}

TMBobj_partialILR_remove_single_obs <- function(sig_obj){
  which_remove <- rowSums(sig_obj$Y == 0) >= (ncol(sig_obj$Y)-1)
  give_subset_samples_TMBobj(sig_obj, rownames(sig_obj$Y)[which_remove])
}

give_subset_samples_TMBobj_2 <- function(TMB_data_obj, selected_rows_obs){
  ## very similar to give_subset_samples_TMBobj, but with samples to keep, not to remove
  TMB_data_obj$Y <- TMB_data_obj$Y[selected_rows_obs,]
  TMB_data_obj$n <- length(selected_rows_obs)
  TMB_data_obj$x <- TMB_data_obj$x[selected_rows_obs,]
  TMB_data_obj$z <- TMB_data_obj$z[selected_rows_obs,]
  TMB_data_obj$z <- TMB_data_obj$z[,colSums(TMB_data_obj$z) >= 1]
  TMB_data_obj$num_individuals <- ncol(TMB_data_obj$z)
  return(TMB_data_obj)
}

make.names_allowing_hyphen <- function(i){
  for(j in unique(i)){
    i[i == j] = paste0(j, '.', 1:sum(i == j))
  }
  i
}

normalise_rw <- function(x){
  if(is.null(dim(x))){
    x/sum(x)
  }else{
    ## normalise row-wise
    sweep(x, 1, rowSums(x), '/')
  }
}

normalise_cl <- function(x){
  if(is.null(dim(x))){
    x/sum(x)
  }else{
    ## normalise col-wise
    t(sweep(x, 2, colSums(x), '/'))
  }
}

give_barplot = function(ct, typedata, simulation=F, title='', legend_on=F, ...){
  require(gridExtra)
  obj = load_PCAWG(ct, typedata, simulation)
  a <- createBarplot(obj$Y[obj$x[,2] == 0,], remove_labels = T, ...)+ggtitle('Early raw')+ guides(shape = guide_legend(override.aes = list(size = 5)))
  b <- createBarplot(obj$Y[obj$x[,2] == 1,], remove_labels = T, ...)+ggtitle('Late raw')+ guides(shape = guide_legend(override.aes = list(size = 5)))
  c <- createBarplot(normalise_rw(obj$Y[obj$x[,2] == 0,]), remove_labels = T, ...)+ggtitle('Early normalised')+ guides(shape = guide_legend(override.aes = list(size = 5)))
  d <- createBarplot(normalise_rw(obj$Y[obj$x[,2] == 1,]), remove_labels = T, ...)+ggtitle('Late normalised')+ guides(shape = guide_legend(override.aes = list(size = 5)))
  
  if(!legend_on){
    a <- a+guides(fill=F)
    b <- b+guides(fill=F)
    c <- c+guides(fill=F)
    d <- d+guides(fill=F)
  }
  grid.arrange(a, b, c, d, top=title)
}

# wrapper_run_TMB_debug = function(object, model = "fullRE_DM", return_report=F, iter.max=150, init_log_lambda = 2, idx_cov_to_fill){
#   dim(object$Y)
#   # sort_columns=T
#   smart_init_vals=T
#   
#   data = object
#   
#   # if(sort_columns){
#   #   data$Y = data$Y[,order(colSums(data$Y), decreasing = F)]
#   # }
#   
#   data$Y = matrix(data$Y, nrow=nrow(data$Y))
#   data$x = (matrix(data$x, ncol=2))
#   
#   d <- ncol(data$Y)
#   n <- ncol(data$z) ## number of INDIVIDUALS, not samples
#   
#   data$num_individuals = n
#   data$lambda_accessory_mat = (cbind(c(rep(1,n),rep(0,n)), c(rep(0,n),rep(1,n))))
#   
#   if(smart_init_vals){
#     require(nnet)
#     .x_multinom = multinom(data$Y ~ data$x[,2])
#     beta_init = t(coef(.x_multinom))
#   }else{
#     beta_init = (matrix(rep(runif(1, min = -4, max = 4), 2*(d-1)),
#                         nrow = 2, byrow=TRUE))
#   }
#   
#   parameters <- list(
#     beta = beta_init,
#     u_large = matrix(runif(min = -0.3, max = 0.3, n = (d-1)*n), nrow=n),
#     logs_sd_RE=rep(1, d-1),
#     # cov_par_RE = rep(1, ((d-1)*(d-1)-(d-1))/2),
#     cov_par_RE = runif(min = -1, max = 1, n = ((d-1)*(d-1)-(d-1))/2),
#     log_lambda = matrix(c(init_log_lambda,init_log_lambda))
#   )
#   if(model == "fullRE_DM"){
#     obj <- MakeADFun(data, parameters, DLL="fullRE_ME_dirichletmultinomial", random = "u_large")
#   }else if(model == "diagRE_DM"){
#     parameters$cov_par_RE = NULL
#     obj <- MakeADFun(data, parameters, DLL="diagRE_ME_dirichletmultinomial", random = "u_large")
#   }else if(model == "fullREDMsinglelambda"){
#     parameters$log_lambda = 2
#     obj <- MakeADFun(data, parameters, DLL="fullRE_dirichletmultinomial_single_lambda", random = "u_large")
#   }else if(model == "diagREDMsinglelambda"){
#     parameters$cov_par_RE = NULL
#     parameters$log_lambda = 2
#     obj <- MakeADFun(data, parameters, DLL="diagRE_dirichletmultinomial_single_lambda", random = "u_large")
#   }else if(model == "sparseRE_DM"){
#     if(is.null(idx_cov_to_fill)){stop("Add <idx_cov_to_fill>")}
#     parameters$cov_par_RE = NULL
#     parameters$cov_RE_part = (rep(1, length(idx_cov_to_fill)))
#     data$idx_params_to_infer = (idx_cov_to_fill)
#     if(length(idx_cov_to_fill) > 1){
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomial", random = "u_large")
#     }else{
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomial_single", random = "u_large")
#     }
#   }else if(model == "sparseRE_DMSL"){
#     if(is.null(idx_cov_to_fill)){stop("Add <idx_cov_to_fill>")}
#     parameters$log_lambda = 2
#     parameters$cov_par_RE = NULL
#     parameters$cov_RE_part = (rep(1, length(idx_cov_to_fill)))
#     data$idx_params_to_infer = (idx_cov_to_fill)
#     if(length(idx_cov_to_fill) > 1){
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomialsinglelambda", random = "u_large")
#     }else{
#       stop()
#     }
#   }else if(model == "sparseRE_DMSL2"){
#     if(is.null(idx_cov_to_fill)){stop("Add <idx_cov_to_fill>")}
#     parameters$log_lambda = 2
#     parameters$cov_par_RE = NULL
#     parameters$cov_RE_part = (rep(1, length(idx_cov_to_fill)))
#     data$idx_params_to_infer = (idx_cov_to_fill)
#     if(length(idx_cov_to_fill) > 1){
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomialsinglelambda2", random = "u_large")
#     }else{
#       stop()
#     }
#   }else if(model == "sparseRE_DMSL2"){
#     if(is.null(idx_cov_to_fill)){stop("Add <idx_cov_to_fill>")}
#     parameters$log_lambda = 2
#     parameters$cov_par_RE = NULL
#     parameters$cov_RE_part = (rep(1, length(idx_cov_to_fill)))
#     data$idx_params_to_infer = (idx_cov_to_fill)
#     if(length(idx_cov_to_fill) > 1){
#       obj <- MakeADFun(data, parameters, DLL="sparseRE_ME_dirichletmultinomialsinglelambda2", random = "u_large")
#     }else{
#       stop()
#     }
#   } else{
#     stop("Incorrect <model>")
#   }
# 
#   opt = nlminb(start = obj$par, obj = obj$fn, gr = obj$gr, iter.max=iter.max, trace=T)
#   
#   if(return_report)  return(sdreport(obj))
#   return(opt)
# }

is_slope = function(v){
  bool_isbetaslope = rep(F, length(v))
  bool_isbetaslope[v == "beta"] = T
  bool_isbetaslope[which(bool_isbetaslope)] = c(F,T)
  bool_isbetaslope
}

give_barplot_from_obj <- function(obj, legend_on=F, legend_bottom=F, nrow_plot=2, only_normalised=F, title=NULL, plot=T, title_facets=NULL, scale_color_manual_vec=NULL, ...){
  if(is.null(title_facets)){
    # title_facets <- c('Early raw', 'Late raw', 'Early normalised', 'Late normalised')
    title_facets <- c('Clonal raw', 'Subclonal raw', 'Clonal normalised', 'Subclonal normalised')
  }
  a <- createBarplot(obj$Y[obj$x[,2] == 0,], remove_labels = T, ...)+ggtitle(title_facets[1])
  b <- createBarplot(obj$Y[obj$x[,2] == 1,], remove_labels = T, ...)+ggtitle(title_facets[2])
  c <- createBarplot(normalise_rw(obj$Y[obj$x[,2] == 0,]), remove_labels = T, ...)+ggtitle(title_facets[3])
  d <- createBarplot(normalise_rw(obj$Y[obj$x[,2] == 1,]), remove_labels = T, ...)+ggtitle(title_facets[4])
  if(!is.null(scale_color_manual_vec)){
    cat('Adding colour scale\n')
    a <- a+scale_fill_manual(values=scale_color_manual_vec)
    b <- b+scale_fill_manual(values=scale_color_manual_vec)
    c <- c+scale_fill_manual(values=scale_color_manual_vec)
    d <- d+scale_fill_manual(values=scale_color_manual_vec)
  }
  if(!legend_on){
    a <- a+guides(fill=F)
    b <- b+guides(fill=F)
    c <- c+guides(fill=F)
    d <- d+guides(fill=F)
  }
  if(legend_bottom){
    a <- a+theme(legend.position='bottom', legend.text = element_text(size = 6), legend.key.size = unit(.2, 'cm'))
    b <- b+theme(legend.position='bottom', legend.text = element_text(size = 6), legend.key.size = unit(.2, 'cm'))
    c <- c+theme(legend.position='bottom', legend.text = element_text(size = 6), legend.key.size = unit(.2, 'cm'))
    d <- d+theme(legend.position='bottom', legend.text = element_text(size = 6), legend.key.size = unit(.2, 'cm'))
  }
  if(only_normalised){
    if(plot){
      grid.arrange(c, d, top=title, nrow=nrow_plot)
    }else{
      cowplot::as_grob(grid.arrange(c, d, top=title, nrow=nrow_plot))
    }
  }else{
    if(plot){
      grid.arrange(a, b, c, d, top=title, nrow=nrow_plot)
    }else{
      cowplot::as_grob(grid.arrange(a, b, c, d, top=title, nrow=nrow_plot))
    }
  }
}

rename_Y_dollar <- function(dataset_obj){
  colnames(dataset_obj$Y) <- gsub(">", "$>$", colnames(dataset_obj$Y))
  dataset_obj
}

give_ranked_plot_simulation = function(tmb_fit_object, data_object, print_plot = T, nreps = 1, model, integer_overdispersion_param){
  
  if(model == 'M'){
    ## theta is always going to be the same. Only replicate the draws from the multinomial
    
    sim_theta = simulate_from_M_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
                                    x_matrix = data_object$x, z_matrix = data_object$z)
    
    sim_counts = t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) )
    
    if(nreps>1){
      sim_counts = replicate(nreps, t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) ), simplify = F)
    }
  }else if(model %in% c('DM', 'DMSL')){
    give_sim_data = function(){
      if(model == 'DM'){
        sim_theta = simulate_from_DM_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
                                         x_matrix = data_object$x, z_matrix = data_object$z, integer_overdispersion_param=integer_overdispersion_param)
      }else if(model == 'DMSL'){
        sim_theta = simulate_from_DMSL_TMB(tmb_fit_object = tmb_fit_object, full_RE = T,
                                           x_matrix = data_object$x, z_matrix = data_object$z, integer_overdispersion_param=integer_overdispersion_param)
      }
      sim_counts = t(sapply(1:nrow(sim_theta), function(i) rmultinom(n = 1, size = sum(data_object$Y[i,]), prob = sim_theta[i,]) ) )
      return(sim_counts)
    }    
    if(nreps>1){
      sim_counts = replicate(nreps, give_sim_data(), simplify = F)
    }else{
      sim_counts = give_sim_data()
    }
  }else{
    stop('Specify a correct model')
  }
  
  stopifnot(all(dim(data_object$Y) == dim(sim_counts)))
  if(print_plot)  plot(sort(data_object$Y), sort(sim_counts))
  
  return(sim_counts)
  
}


simulate_from_M_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = softmax(cbind(logRmat, 0))
  }else{
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
  }
  return(sim_thetas)
}

simulate_from_DM_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, integer_overdispersion_param){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  overdispersion_lambda = integer_overdispersion_param*exp(python_like_select_name(tmb_fit_object$par.fixed, "log_lambda")[x_matrix[,2]+1])
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* softmax(c(logRmat[l,], 0)))))
  }else{
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* sim_thetas[l,])))
  }
  return(sim_thetas)
}

simulate_from_DMSL_TMB = function(tmb_fit_object, full_RE=T, x_matrix, z_matrix, integer_overdispersion_param){
  dmin1 = length(python_like_select_name(tmb_fit_object$par.fixed, 'beta'))/2
  overdispersion_lambda = rep(integer_overdispersion_param*exp(python_like_select_name(tmb_fit_object$par.fixed, "log_lambda")), nrow(x_matrix))
  if(full_RE){
    re_mat = re_vector_to_matrix(tmb_fit_object$par.random, dmin1)
    ntimes2 = nrow(z_matrix)
    logRmat = z_matrix %*% re_mat + 
      x_matrix %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2)
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* softmax(c(logRmat[l,], 0)))))
  }else{
    ## univariate RE
    sim_thetas = softmax(cbind(sapply(1:dmin1,
                                      function(some_dummy_idx) give_z_matrix(length(tmb_fit_object$par.random) * 2) %*% tmb_fit_object$par.random) +
                                 give_x_matrix(length(tmb_fit_object$par.random) * 2) %*% matrix(python_like_select_name(tmb_fit_object$par.fixed, 'beta'), nrow=2), 0))
    sim_thetas = t(sapply(1:nrow(logRmat), function(l) MCMCpack::rdirichlet(1, overdispersion_lambda[l]* sim_thetas[l,])))
  }
  return(sim_thetas)
}

give_interval_plots = function(df_rank, data_object, loglog=F){
  xx = melt(df_rank, id.vars=c('sorted_value', 'rank_number'))
  xx_summary = xx %>% group_by(rank_number) %>% mutate(min_interval=quantile(sorted_value, probs = c(0.025)),
                                                       max_interval=quantile(sorted_value, probs = c(0.975)),
                                                       mean=mean(sorted_value))
  xx_summary = xx_summary[!(duplicated(xx_summary[,c('rank_number', 'min_interval', 'max_interval')])),c('rank_number', 'min_interval', 'max_interval', 'mean')]
  a = ggplot(cbind.data.frame(xx_summary, sorted_true=sort(data_object$Y)), aes(x=sorted_true, ymin= min_interval, ymax=max_interval))+
    geom_abline(slope = 1, intercept = 0, lty='dashed')+
    geom_ribbon(fill = "red", alpha=0.2)+
    geom_point(aes(x=sorted_true, y=mean), size=0.4)+
    geom_line(aes(x=sorted_true, y=mean))+
    labs(x='Observed ranked value', y='Mean simulated ranked value')
  if(loglog){
    a = a+ scale_y_continuous(trans = "log")+scale_x_continuous(trans = "log")
  }
  
  return(a)
}

give_interval_plots_2 = function(df_rank, data_object,loglog=F, title, theme_bw=T){
  xx = melt(df_rank, id.vars=c('sorted_value', 'rank_number'))
  xx_summary = xx %>% group_by(rank_number) %>% mutate(min_interval=quantile(sorted_value, probs = c(0.025)),
                                                       max_interval=quantile(sorted_value, probs = c(0.975)),
                                                       mean=mean(sorted_value))
  xx_summary = xx_summary[!(duplicated(xx_summary[,c('rank_number', 'min_interval', 'max_interval')])),c('rank_number', 'min_interval', 'max_interval', 'mean')]
  xx_summary <- cbind.data.frame(xx_summary, sorted_true=as.vector(data_object$Y))
  xx_summary[,'rank_true'] = order(xx_summary$sorted_true)
  xx_summary$col = apply(xx_summary, 1, function(i) (i['sorted_true'] > i['min_interval']) & (i['sorted_true'] < i['max_interval']) )
  a  = ggplot(xx_summary, aes(x=sorted_true,
                              ymin= min_interval, ymax=max_interval))+
    geom_abline(slope = 1, intercept = 0, lty='dashed')+
    geom_ribbon(fill = "red", alpha=0.2)+
    geom_point(aes(x=sorted_true, y=mean, col=col), size=0.4)+
    geom_line(aes(x=sorted_true, y=mean))+
    labs(x='Observed value', y='Mean simulated value')+
    ggtitle(title, subtitle=paste0(paste0(names(table(xx_summary$col)), ':', table(xx_summary$col)), collapse ='; '))+
    theme(legend.position = "bottom")
  
  if(theme_bw){
    a = a+theme_bw()+theme(legend.position = "bottom")
  }
  
  if(loglog){a=a+ scale_y_continuous(trans = "log")+scale_x_continuous(trans = "log") }
  
  return(a)
}

give_betas <- function(TMB_obj){
  matrix(python_like_select_name(TMB_obj$par.fixed, 'beta'), nrow=2)
}



plot_betas <- function(TMB_obj, names_cats=NULL, rotate_axis=T, theme_bw=T, remove_SBS=T, only_slope=F, return_df=F, plot=T,
                       line_zero=T, add_confint=F, return_plot=T, return_ggplot=F, title=NULL, add_median=F, sort_by_slope=F,
                       size_title=NULL, size_logR=NULL, xlab=NULL, betas_are_only_slope=F, input_is_summary=F){
  if(typeof(TMB_obj) == 'character'){
    .summary_betas <- NA
    if(theme_bw){
      plt <- ggplot()+theme_bw()
      if(plot) print(plt)
    }else{
      plt <- ggplot()
      if(plot) print(plt)
    }
  }else{
    if(!input_is_summary){
      .summary_betas <- summary(TMB_obj)
    }else{
      .summary_betas <- TMB_obj
    }
    if(betas_are_only_slope){
      .extra_beta <- matrix(NA, nrow=nrow(python_like_select_rownames(.summary_betas, 'beta')), ncol=2)
      .summary_betas
      .extra_beta <- do.call('rbind', lapply(1:nrow(.extra_beta), function(i) rbind(.extra_beta[i,], .summary_betas[i,])))
      rownames(.extra_beta) <- rep('beta', nrow(.extra_beta))
      .summary_betas <- .summary_betas[!rownames(.summary_betas) == 'beta',]
      .summary_betas <- rbind(.extra_beta, .summary_betas)
    }
    .summary_betas <- cbind.data.frame(python_like_select_rownames(.summary_betas, 'beta'),
                                       type_beta=rep(c('Intercept', 'Slope')),
                                       LogR=rep(1:(nrow(python_like_select_rownames(.summary_betas, 'beta'))/2), each=2))
    if(only_slope){
      .summary_betas <- .summary_betas[.summary_betas$type_beta == 'Slope',]
    }
    
    if(!is.null(names_cats)){
      if(remove_SBS){
        names_cats <- gsub("SBS", "", names_cats) 
      }
      if(length(unique(.summary_betas$LogR)) != length(names_cats)){
        stop('Number of beta slope/intercept pairs should be the same as the length of the name of the categories')
      }
      .summary_betas$LogR = names_cats[.summary_betas$LogR]
    }
    if(sort_by_slope){
      .summary_betas$LogR <- factor( .summary_betas$LogR,
                                     levels=.summary_betas[.summary_betas$type_beta == 'Slope','LogR'][order(.summary_betas[.summary_betas$type_beta == 'Slope','Estimate'])])
    }
    plt <- ggplot(.summary_betas, aes(x=LogR, y=`Estimate`))
    
    if(line_zero) plt <- plt + geom_hline(yintercept = 0, lty='dashed', col='blue')
    if(add_median) { plt <- plt +
      geom_hline(yintercept = median(c(0,.summary_betas$Estimate[.summary_betas$type_beta == 'Slope'])),
                 lty='dashed', col='red') }
    
    plt <- plt +
      geom_point()+
      geom_errorbar(aes(ymin=`Estimate`-`Std. Error`, ymax=`Estimate`+`Std. Error`), width=.1)+
      facet_wrap(.~type_beta, scales = "free") #ggtitle('Slopes')
    
    if(theme_bw){
      plt <- plt + theme_bw()
    }
    
    if(add_confint){
      confints <- cbind(.summary_betas, confint=t(give_confidence_interval(.summary_betas[,'Estimate'], .summary_betas[,'Std. Error'])))
      plt <- plt+
        geom_errorbar(data = confints, aes(ymin=confint.1, ymax=confint.2), width=.1 ,col='blue', alpha=0.6)
      
    }
    
    if(rotate_axis){
      plt <- plt + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
    if(!is.null(size_title)){
      plt <- plt + theme(plot.title = element_text(size=size_title))
    }
    
    if(!is.null(size_logR)){
      plt <- plt + theme(axis.text.x = element_text(size=size_logR))
    }
    
    if(!is.null(title)){
      plt <- plt + ggtitle(title)
    }
    
    if(!input_is_summary){
      if(!TMB_obj$pdHess){
        plt <- plt + annotate("text", x = -Inf, y=Inf, label="not PD", vjust=1, hjust=-.2)+geom_point(col='red')
      }
    }
    if(plot) print(plt)
      
  }
  
  if(!is.null(xlab)){
    cat('Changing x axis label\n')
    plt <- plt+labs(x=xlab)
  }
  
  if(return_df){
    .summary_betas
  }else{
    if(return_plot & return_df){stop('<return_plot=T> and <return_df=T> are incompatible')}
    plot_list <- list(plt)
    class(plot_list) <- c("quiet_list", class(plot_list))
    if(return_plot){
      return(cowplot::as_grob(plt))
    }else if(return_ggplot){
      return(plt)
    }
  }
}

compare_betas_tmb <- function(tmb_obj_1, names_cats1, tmb_obj_2, names_cats2, names_groups=c('Group 1', 'Group 2'),
                              include_missing_as_inf=F){
  betas1 <- softmax(cbind(matrix(python_like_select_name(tmb_obj_1$par.fixed, 'beta'),  nrow=2), 0))
  betas2 <- softmax(cbind(matrix(python_like_select_name(tmb_obj_2$par.fixed, 'beta'), nrow=2), 0))
  betas1
  colnames(betas1) <- names_cats1
  colnames(betas2) <- names_cats2
  
  if(include_missing_as_inf){
    betas2 <- betas2[,(match(colnames(betas1), colnames(betas2)))]
    colnames(betas2) <- colnames(betas1)
    betas2[is.na(betas2)] <- Inf
    betas1 <- betas1[,match(colnames(betas2), colnames(betas1))]    
  }else{
    betas2 <- betas2[,remove_na(match(colnames(betas1), colnames(betas2)))]
    betas1 <- betas1[,match(colnames(betas2), colnames(betas1))]
  }
  
  .x <- cbind.data.frame(obj1=melt(betas1),
                         obj2=melt(betas2))
  .x$obj1.Var1 <- c('Intercept', 'Slope')[.x$obj1.Var1]
  .x$col <- (is.infinite(.x$obj1.value) | is.infinite(.x$obj2.value))
  ggplot(.x, aes(x=obj1.value, y=obj2.value, label=gsub('SBS', '', obj1.Var2),
                 alpha=1-col))+
    geom_abline(slope = 1, intercept = 0, lty='dashed')+
    geom_point()+
    facet_wrap(.~obj1.Var1, scales = "free")+theme_bw()+geom_label_repel()+
    labs(x=paste0('Softmax beta ', names_groups[1]),
         y=paste0('Softmax beta ', names_groups[2]))
}

remove_na <- function(i) i[!is.na(i)]

nan_to_zero <- function(m){
  m[is.nan(m)] <- 0
}

rm_na_rows <- function(m){
  m[!(apply(m, 1, function(i) any(is.na(i)))),]
}

give_ternary <- function(probs, add_par=T, opacity=0.2, col='red', pch='.', cex=5, ...){
  require(Ternary)
  mar=c(0,0,0,0)
  TernaryPlot(atip = colnames(probs)[1], btip = colnames(probs)[2], ctip = colnames(probs)[3], grid.lines = 0, grid.col = NULL, ...)
  dens <- TernaryDensity(probs, resolution = 10L)
  
  cls_legend = rbind(viridisLite::viridis(48L, alpha = 0.6),
                     seq(from = 0, to = 47, by=1))
  TernaryPoints(probs, col = alpha(col, opacity), pch = pch, cex=cex)
}

give_ternary_v2 <- function(probs, add_par=T, opacity=0.2, col='red', pch='.', cex=5, legend_off=F, ...){
  par(mar=c(0,0,0,0))
  TernaryPlot(atip = colnames(probs)[1], btip = colnames(probs)[2], ctip = colnames(probs)[3],
              grid.lines = 0, grid.col = NULL)
  dens <- TernaryDensity(probs, resolution = 10L)
  
  cls_legend = rbind(viridisLite::viridis(48L, alpha = 0.6),
                     seq(from = 0, to = 47, by=1))
  if(!legend_off){
    legend(x=-0.6,y=1.08,
           fill = cls_legend[1,][c(T,F,F,F,F)],
           legend = round(as.numeric(cls_legend[2,][c(T,F,F,F,F)])/sum(dens['z',]), 2), ncol=5,
           y.intersp=0.8,x.intersp=0.5,text.width=0.1, cex=0.9, bty = "n")
  }
  ColourTernary(dens)
  TernaryPoints(probs, col = 'red', pch = '.', cex=5)
  TernaryDensityContour(probs, resolution = 30L)
}

object_to_cov_heatmap <- function(tmb_obj, d, names_cov=NULL){
  cov <- L_to_cov(python_like_select_name(tmb_obj$par.fixed, 'cov_RE'),
                  d = d)
  if(is.null(names_cov)){
    names_cov <- paste0('ALR', 1:d)
  }
  colnames(cov) <- rownames(cov) <- names_cov
  pheatmap(cov)
}


give_amalgamation <- function(i, list_amalgamations){
  new_mat = sapply(list_amalgamations, function(j){
    grouped_exp <- i[,colnames(i) %in% j]
    if(!is.null(ncol(grouped_exp))){
      rowSums(grouped_exp)
    }else{
      grouped_exp
    }
  })
  new_mat
}

plot_lambdas <- function(TMB_obj, return_df=F, plot=T){
  
  lambdas_df <- python_like_select_rownames(summary(TMB_obj), 'log_lambda')
  if(!is.null(dim(lambdas_df))){
    .summary_lambda <- cbind.data.frame(data.frame(lambdas_df),
                                        name=c('Lambda 1', 'Lambda 2'))
  }else{
    .summary_lambda <- cbind.data.frame(t(lambdas_df), name='Lambda 1')
    colnames(.summary_lambda) <- make.names(colnames(.summary_lambda))
  }
  
  plt <- (ggplot(.summary_lambda, aes(x=name, y=`Estimate`))+
            geom_point()+
            geom_errorbar(aes(ymin=`Estimate`-`Std..Error`, ymax=`Estimate`+`Std..Error`), width=.1)+theme_bw())
  if(plot){
    print(plt)
  }
  
  if(return_df){
    return(.summary_lambda)
  }else{
    return(plt)
  }
}

plot_estimates_TMB <- function(TMB_obj, parameter_name, return_df=F, plot=T, verbatim=T){
  if(verbatim) cat('Consider the other functions <plot_betas> and <plot_lambdas>\n')
  
  
  parameter_df <- python_like_select_rownames(summary(TMB_obj), parameter_name)
  if(length(parameter_df) == 0){
    ## there is no parameter
    return(cbind.data.frame(Estimate=NA, `Std..Error`=NA, name=NA))
  }else{
    if(!is.null(dim(parameter_df))){
      .summary_param <- cbind.data.frame(data.frame(parameter_df),
                                         name=paste0(parameter_name, 1:nrow(parameter_df)))
    }else{
      ## single parameter
      .summary_param <- cbind.data.frame(t(parameter_df), name=paste0(parameter_name, '1'))
      colnames(.summary_param) <- make.names(colnames(.summary_param))
    }
    
    plt <- (ggplot(.summary_param, aes(x=name, y=`Estimate`))+
              geom_point()+
              geom_errorbar(aes(ymin=`Estimate`-`Std..Error`, ymax=`Estimate`+`Std..Error`), width=.1)+theme_bw())
    if(plot){
      print(plt)
    }
    
    if(return_df){
      return(.summary_param)
    }else{
      return(plt)
    }
  }
}

give_length_cov <- function(dim1_covmat){
  ((dim1_covmat**2) - dim1_covmat )/2
}


give_sim_from_estimates <- function(ct, typedata = "signatures", sigs_to_remove="", model="sparseRE_DM",
                                    bool_nonexo=TRUE, bool_give_PCA, sig_of_interest='SBS8',
                                    path_to_data= "../../../data/", tmb_object=NULL, obj_data=NULL,
                                    nrow_pca_plot=2, integer_overdispersion_param=1){
  warning('<fullRE_DMSL> and <fullRE_DM> used be to used interchangeably')
  
  if(is.null(tmb_object)){
    if(model == "fullRE_M"){
      if(!bool_nonexo)    list_estimates <- fullRE_M
      if(bool_nonexo)    list_estimates <- fullRE_M_nonexo
    }else if(model == "fullRE_DMSL"){
      if(!bool_nonexo)    list_estimates <- fullRE_DMSL
      if(bool_nonexo)    list_estimates <- fullRE_DMSL_nonexo
    }else if(model == "sparseRE_DM"){
      if(bool_nonexo)    list_estimates <- sparseRE_DMSL_nonexo
    }else{
      stop('Add <tmb_object> or specify another model')
    }
  }else{
    list_estimates <- list()
    list_estimates[[ct]] <- tmb_object
  }
  
  if(is.null(obj_data)){
    warning('WARNING! Here I am sorting the columns of TMB, I should not in all cases. Specify <obj_data> if needed')
    obj_data <- sort_columns_TMB(give_subset_sigs_TMBobj(load_PCAWG(ct = ct, typedata = typedata, path_to_data =path_to_data),
                                                         sigs_to_remove = sigs_to_remove))
  }
  dmin1 <- ncol(obj_data$Y)-1
  
  if(model %in% c("fullRE_M", "fullRE_DM", "fullRE_DMSL")){
    ### implementing var_vec as though they were sd ###
    cov_mat = L_to_cov(python_like_select_name(list_estimates[[ct]]$par.fixed, 'cov_par_RE'), d=dmin1)
  }else if(model == "sparseRE_DM"){
    cov_vec = rep(0, (dmin1**2-dmin1)/2)
    cov_vec[as.numeric(strsplit(subset_sigs_sparse_cov_idx_nonexo[subset_sigs_sparse_cov_idx_nonexo$V1 == ct,"V2"], ',')[[1]])] = python_like_select_name(list_estimates[[ct]]$par.fixed, 'cov_RE_part')
    cov_mat = L_to_cov(cov_vec, d=dmin1)
  }else if(model %in% c("diagRE_M", "diagRE_DMSL", "diagRE_DMDL")){
    cov_vec = rep(0, give_length_cov(length(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))))
    cov_mat = L_to_cov(cov_vec, d=dmin1)
  }else{
    stop('Specify correct <model>')
  }
  var_vec = exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'logs_sd_RE'))**2
  diag(cov_mat) <- var_vec
  
  beta_mat = matrix(python_like_select_name(list_estimates[[ct]]$par.fixed, 'beta'), nrow=2)
  
  n_sim = 1000
  x_sim = cbind(1, rep(c(0,1), n_sim))
  u_sim = mvtnorm::rmvnorm(n = n_sim, mean = rep(0,dmin1), sigma = cov_mat)
  
  theta = x_sim %*% beta_mat + (give_z_matrix(n_sim*2)) %*% u_sim
  
  if(model %in% c('sparseRE_DM', 'fullRE_DM', 'fullRE_DMSL', 'diagRE_DMSL')){
    ## single overdispersion parameter
    alpha = softmax(cbind(theta, 0))*exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'log_lambda'))*integer_overdispersion_param
  }else if(model %in% c("fullRE_M")){
    ## no overdispersion parameter
    alpha = softmax(cbind(theta, 0))
  }else if(model %in% c('fullRE_DMDL', 'diagRE_DMDL')){
    ## two overdispersion parameters
    alpha = softmax(cbind(theta, 0))*c(exp(python_like_select_name(list_estimates[[ct]]$par.fixed, 'log_lambda'))[rep(c(1,2), each=n_sim)]*integer_overdispersion_param)
  }else{
    stop('Check softmax step')
  }
  
  if(model %in% c('sparseRE_DM', 'fullRE_DM', 'fullRE_DMSL', 'diagRE_DMSL', 'diagRE_DMDL')){
    probs = t(apply(alpha, 1, MCMCpack::rdirichlet, n=1))
  }else if(model %in% c("fullRE_M")){
    probs = alpha
  }else{
    stop('Check probabilities step')
  }
  
  probs_obs = normalise_rw(obj_data$Y)
  all_probs = rbind(probs_obs, probs)
  
  if(bool_give_PCA){
    # pca <- prcomp(as(compositions::clr(all_probs), 'matrix'))
    pca <- prcomp(all_probs)
    
    if(! (sig_of_interest %in% colnames(all_probs))){stop('Specify a <sig_of_interest> present in the dataset')}
    
    df_pca <- cbind.data.frame(pca=pca$x[,1:2], col=c(rep('Observed', nrow(probs_obs)), rep('Simulated',nrow(probs))),
                               sig_of_interest=all_probs[,sig_of_interest],
                               group=c(rep(c('early','late')[obj_data$x[,2]+1]),
                                       rep(c('early','late'), n_sim)))
    return(list(df_pca, ggplot(df_pca, aes(x=pca.PC1, y=pca.PC2, col=sig_of_interest))+labs(col=sig_of_interest)+
                  geom_point(alpha=0.7)+facet_wrap(.~interaction(col,group),nrow=nrow_pca_plot)+theme_bw()))
  }else{
    return(all_probs)
  }
  
}

select_self <- function(i) i[i]

vector_cats_to_logR <- function(i){paste0(i[-length(i)], '/', i[length(i)])}

non_duplicated_rows <- function(i){
  rownames(i)[duplicated(rownames(i))] = paste0(rownames(i)[duplicated(rownames(i))], "_2")
  i
}

give_all_correlations <- function(j){
  outer(1:nrow(j), 1:nrow(j), Vectorize(function(i,k){cor(j[i,], j[k,])}))
}

give_perturbation_TMBobj = function(exposures_cancertype_obj, addone=F){
  ## adapted from the function give_perturbation_ROOSigs_alt
  
  if(addone){
    exposures_cancertype_obj$Y = exposures_cancertype_obj$Y + 1
  }
  
  ## for each individual, compute the perturbation
  apply(exposures_cancertype_obj$z, 2, function(i){
    ## if there are two individuals
    if(sum(i) == 2){
      exposures_cancertype_obj$Y[which(i == 1)[1],] / exposures_cancertype_obj$Y[which(i == 1)[2],]
    }
  })
  
}

give_totalperturbation_TMBobj = function(exposures_cancertype_obj, addone=F){
  ## adapted from the function give_total_perturbation
  
  pert <- give_perturbation_TMBobj(exposures_cancertype_obj, addone=addone)
  
  mean(apply(pert, 2, function(i) sqrt(sum((i-1/(ncol(exposures_cancertype_obj$Y)))**2))))
  
}

give_totalperturbation_TMBobj_sigaverage = function(exposures_cancertype_obj, addone=F){
  ## adapted from the function give_total_perturbation
  
  pert <- give_perturbation_TMBobj(exposures_cancertype_obj, addone=addone)
  
  mean(apply(pert, 1, function(i) sqrt(sum((i-1/(ncol(exposures_cancertype_obj$Y)))**2))))
  
}

give_weightedtotalperturbation_TMBobj = function(exposures_cancertype_obj, addone=F){
  
  pert <- give_perturbation_TMBobj(exposures_cancertype_obj, addone=addone)
  
  perts <- apply(pert, 2, function(i) sqrt(sum((i-1/(ncol(exposures_cancertype_obj$Y)))**2)))
  weights <- as.vector(apply(exposures_cancertype_obj$z, 2, function(i){
    ## if there are two individuals
    if(sum(i) == 2){
      log(sum(rowSums(exposures_cancertype_obj$Y[which(i == 1),])))
    }
  }))
  if(length(perts) != length(weights)){
    stop('Weights and perturbations do not have the same length')
  }
  weights <- weights/sum(weights) # normalise
  sum(perts * weights)/nrow(exposures_cancertype_obj$Y)
  
}

wrapper_run_HMP_Xdc.sevsample <- function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  return(HMP::Xdc.sevsample(x)$`p value`)
}


wrapper_run_HMP_Xmcupo.sevsample <- function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  return(HMP::Xmcupo.sevsample(x)$`p value`)
}

plot_ternary <- function(x, legend_on=T, plot_points=T, ...){
  require(Ternary)
  
  if(ncol(x) != 3){stop('Number of columns must be three. Create a subcomposition or amalgamation if needed.')}
  TernaryPlot(atip = colnames(x)[1], btip = colnames(x)[2], ctip = colnames(x)[3],
              grid.lines = 0, grid.col = NULL, ...)
  dens <- TernaryDensity(x, resolution = 10L)
  
  cls_legend = rbind(viridisLite::viridis(48L, alpha = 0.6),
                     seq(from = 0, to = 47, by=1))
  if(legend_on){
    legend(x=-0.4,y=1.08,
           fill = cls_legend[1,][c(T,F,F,F,F)],
           legend = round(as.numeric(cls_legend[2,][c(T,F,F,F,F)])/sum(dens['z',]), 2), ncol=5,
           y.intersp=0.8,x.intersp=0.5,text.width=0.1, cex=0.9, bty = "n")
  }
  ColourTernary(dens)
  if(plot_points)  TernaryPoints(x, col = 'red', pch = '.', cex=5)
  TernaryDensityContour(x, resolution = 30L)
}

split_matrix_in_half <- function(x){
  list(x[1:(nrow(x)/2),],
       x[(1+(nrow(x)/2)):nrow(x),])
}

createBarplot <- function(matrix_exposures, angle_rotation_axis = 0, order_labels=NULL,
                          remove_labels=FALSE, levels_signatures=NULL, includeMelt=NULL,
                          Melt=NULL, verbose=TRUE, arg_title='Signature', reorder_sigs=NULL){
  #' error due to it not being a matrix:  
  #'    No id variables; using all as measure variables
  #'    Rerun with Debug
  #'    Error in `[.data.frame`(.mat, , "Var1") : undefined columns selected \\
  #' (use tomatrix())
  
  
  if(is.null(colnames(matrix_exposures))) stop('columns must have names')
  require(reshape2)
  require(ggplot2)
  library(RColorBrewer)
  if(verbose){
    cat(paste0('Creating plot... it might take some time if the data are large. Number of samples: ', nrow(matrix_exposures), '\n'))
  }
  
  if( (!is.null(order_labels) & typeof(order_labels) == "logical")){if(!order_labels){cat('WARNING: Order labels is either a vector with desired order or NULL, not bool')}}
  
  if(!is.null(levels_signatures)){
    library(ggthemes)
    ggthemes_data$economist
  }else{
    levels_signatures <- colnames(matrix_exposures)
  }
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  col_vector <- c(col_vector[c(T,F)], col_vector[c(F,T)])
  
  myColors <- col_vector[1:length(levels_signatures)]
  names(myColors) <- levels_signatures
  myColors <- myColors[levels_signatures %in% unique(colnames(matrix_exposures))]
  if(is.null(order_labels)) order_labels = rownames(matrix_exposures)
  if(!is.null(includeMelt)){
    cat('For whatever reason sometimes the melt does not work. Here it is passed as argument.')
    .mat <- Melt
  }else{
    .mat <- melt(matrix_exposures)
  }
  .mat[,'Var1'] <- factor(.mat[,'Var1'], levels=order_labels)
  .mat[,'Var2'] <- factor(.mat[,'Var2'], levels=levels_signatures)
  ###rownames(.mat) <- rownames(matrix_exposures) ### new
  
  if(!is.null(reorder_sigs)){
    ### if we want signatures to be re-ordered, but for them to have a common colour scheme,
    ## use reorder_sigs together with levels_signatures, where levels_signatures creates a 
    ## common colour scheme, and reorder_sigs changes the order of the bars
    levels_signatures <- reorder_sigs
  }
  
  if(!remove_labels){
    if(!is.null(levels_signatures)){
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures[levels_signatures %in% unique(colnames(matrix_exposures))])))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        #theme(axis.title.x=element_blank(),
        #      axis.text.x=element_blank(),
        #      axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title=arg_title))+
        scale_fill_manual(name = "grp",values = myColors)
    }else{
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        #theme(axis.title.x=element_blank(),
        #      axis.text.x=element_blank(),
        #      axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title=arg_title))
    }
  }else{
    if(!is.null(levels_signatures)){
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title=arg_title))+
        scale_fill_manual(name = "grp",values = myColors)
    }else{
      ggplot(.mat, aes(x=Var1, y=value, fill=factor(Var2, levels=levels_signatures)))+
        geom_bar(stat = 'identity')+
        theme(axis.text.x = element_text(angle = angle_rotation_axis, hjust = 1))+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title=arg_title))
    }
  }
}

give_min_pert <- function(idx_sp, list_runs=diagRE_DMDL_nonexo_SP, logR_names_vec=logR_nonexo_notsorted_SP,
                          df_betas=NULL){
  require(dplyr)
  
  if(!is.null(df_betas)){
    ## get betas from df
    .summary_betas_slope_SP <- df_betas
    .slopes_minpert_SP <- .summary_betas_slope_SP[,1]
  }else{
    ## get betas from TMB output
    .betas_SP <- data.frame(plot_betas(list_runs[[idx_sp]], names_cats= logR_names_vec[[idx_sp]],
                                       return_df=T, plot=F))
    
    .slopes_minpert_SP <- .betas_SP %>% dplyr::filter(type_beta == "Slope") %>% dplyr::select(Estimate) %>% unlist()
    # print(.slopes_minpert_SP)
    ## check if the CI of the betas touches this median value
    .summary_betas_slope_SP <- python_like_select_rownames(summary(list_runs[[idx_sp]]), 'beta')[c(F,T),]
    nrow(.summary_betas_slope_SP)
  }
  
  minimal_change_baseline <- median(c(.slopes_minpert_SP, 0))
  # print(.summary_betas_slope_SP)
  # print(logR_nonexo_notsorted_SP[[idx_sp]])
  # print(dim(.summary_betas_slope_SP))
  if(!is.null(dim(.summary_betas_slope_SP))){
    .params_in_ci <- give_params_in_CI(vec_est=.summary_betas_slope_SP[,1],
                                       vec_stderr=.summary_betas_slope_SP[,2],
                                       vec_true=rep(minimal_change_baseline, nrow(.summary_betas_slope_SP)))
  }else{
    .params_in_ci <- give_params_in_CI(vec_est=.summary_betas_slope_SP[1],
                                       vec_stderr=.summary_betas_slope_SP[2],
                                       vec_true=minimal_change_baseline)
  }
  .params_in_ci <- sapply(1:length(.params_in_ci), function(i){
    ## for the ones in which there is a change, say whether it's up- or down-regulated
    if(!.params_in_ci[i]){
      ## if there is a change: not in confidence interval
      if(is.null(dim(.summary_betas_slope_SP))){
        ## one-dim
        if(.summary_betas_slope_SP[1] > minimal_change_baseline){
          'increase'
        }else{
          'decrease'
        }
      }else{
        ## multi-dim
        if(.summary_betas_slope_SP[i,1] > minimal_change_baseline){
          'increase'
        }else{
          'decrease'
        }
      }
    }else{
      'FALSE'
    }
  })
  
  if(!is.null(df_betas)){
    .baseline <- NA
  }else{
    names(.params_in_ci) <- sapply(logR_names_vec[[idx_sp]], function(i) strsplit(i, '/')[[1]][1])
    .baseline <- strsplit(logR_names_vec[[idx_sp]][[1]], '/')[[1]][2]
  }
  return(list(betas_perturbed=.params_in_ci, baseline=.baseline))
}

comparison_betas_models <- function(model_fullRE_DMSL_list, model_diagRE_DMSL_list, model_fullRE_M_list ){
  
  .x <- do.call('rbind.data.frame', lapply(enough_samples, function(ct){
    x_beta_fullRE_DMSL <- try(python_like_select_name(model_fullRE_DMSL_list[[ct]]$par.fixed, "beta"))
    x_beta_diagRE_DMSL <- try(python_like_select_name(model_diagRE_DMSL_list[[ct]]$par.fixed, "beta"))
    x_beta_fullRE_M <- try(python_like_select_name(model_fullRE_M_list[[ct]]$par.fixed, "beta"))
    if( (length(x_beta_fullRE_DMSL) != length(x_beta_diagRE_DMSL)) | (length(x_beta_fullRE_DMSL) != length(x_beta_fullRE_M)) ){
      ## if we don't have results for any, remove from the analysis
      list_betas <- list(x_beta_fullRE_DMSL, x_beta_diagRE_DMSL, x_beta_fullRE_M)
      typeofs_of_betas <- sapply(list_betas, typeof)
      if( all(typeofs_of_betas == "character")  ){
        return(NULL)
      }else{
        ## if we do have results for some, replace the error message by an NA string
        ## replace using the length of the first double entry
        
        if(typeofs_of_betas[1] == "character"){
          x_beta_fullRE_DMSL <- rep(NA, length(list_betas[[which(typeofs_of_betas == "double")[1]]]))
        }
        
        if(typeofs_of_betas[2] == "character"){
          x_beta_diagRE_DMSL <- rep(NA, length(list_betas[[which(typeofs_of_betas == "double")[1]]]))
        }
        
        if(typeofs_of_betas[3] == "character"){
          x_beta_fullRE_M <- rep(NA, length(list_betas[[which(typeofs_of_betas == "double")[1]]]))
        }
        
        if( (length(x_beta_fullRE_DMSL) != length(x_beta_diagRE_DMSL)) | (length(x_beta_fullRE_DMSL) != length(x_beta_fullRE_M)) ){
          warning(paste0(ct, ': the number of log-ratios is not consistent'))
          return(NULL)
        }
      }
      
    }
    cbind.data.frame(rbind.data.frame(
      cbind.data.frame(fullRE_DMSL=select_slope_2(x_beta_fullRE_DMSL),
                       diagRE_DMSL=select_slope_2(x_beta_diagRE_DMSL),
                       fullRE_M=select_slope_2(x_beta_fullRE_M),
                       beta_type='Slope'),
      cbind.data.frame(fullRE_DMSL=select_intercept(x_beta_fullRE_DMSL),
                       diagRE_DMSL=select_intercept(x_beta_diagRE_DMSL),
                       fullRE_M=select_intercept(x_beta_fullRE_M),
                       beta_type='Intercept')),
      ct=ct)
  }
  ))
  
  ## if something hasn't converged, remove the value
  bad_ct_fullRE_DMSL <- c(enough_samples[sapply(enough_samples,
                                                function(ct) (typeof(model_fullRE_DMSL_list[[ct]])) == 'character')],
                          enough_samples[sapply(enough_samples,
                                                function(ct)  try(model_fullRE_DMSL_list[[ct]]$pdHess)) != "TRUE"])
  bad_ct_diagRE_DMSL <- c(enough_samples[sapply(enough_samples,
                                                function(ct) (typeof(model_diagRE_DMSL_list[[ct]])) == 'character')],
                          enough_samples[sapply(enough_samples, function(ct)  try(model_diagRE_DMSL_list[[ct]]$pdHess)) != "TRUE"])
  bad_ct_fullRE_M <- c(enough_samples[sapply(enough_samples,
                                             function(ct) (typeof(model_fullRE_M_list[[ct]])) == 'character')],
                       enough_samples[sapply(enough_samples, function(ct)  try(model_fullRE_M_list[[ct]]$pdHess)) != "TRUE"])
  
  .x$fullRE_DMSL[(.x$ct %in% bad_ct_fullRE_DMSL)] = NA
  .x$diagRE_DMSL[(.x$ct %in% bad_ct_diagRE_DMSL)] = NA
  .x$fullRE_M[(.x$ct %in% bad_ct_fullRE_M)] = NA
  
  .x$ct2=renaming_pcawg[,2][match(.x$ct, renaming_pcawg[,1])]
  
  return(.x)
  
}

comparison_betas_models2 <- function(model_fullRE_DMSL_list, model_diagRE_DMDL_list, model_fullRE_M_list ){
  
  ## with model_diagRE_DMDL_list instead of model_diagRE_DMSL_list
  
  .x <- do.call('rbind.data.frame', lapply(enough_samples, function(ct){
    x_beta_fullRE_DMSL <- try(python_like_select_name(model_fullRE_DMSL_list[[ct]]$par.fixed, "beta"))
    x_beta_diagRE_DMDL <- try(python_like_select_name(model_diagRE_DMDL_list[[ct]]$par.fixed, "beta"))
    x_beta_fullRE_M <- try(python_like_select_name(model_fullRE_M_list[[ct]]$par.fixed, "beta"))
    if( (length(x_beta_fullRE_DMSL) != length(x_beta_diagRE_DMDL)) | (length(x_beta_fullRE_DMSL) != length(x_beta_fullRE_M)) ){
      ## if we don't have results for any, remove from the analysis
      list_betas <- list(x_beta_fullRE_DMSL, x_beta_diagRE_DMDL, x_beta_fullRE_M)
      typeofs_of_betas <- sapply(list_betas, typeof)
      if( all(typeofs_of_betas == "character")  ){
        return(NULL)
      }else{
        ## if we do have results for some, replace the error message by an NA string
        ## replace using the length of the first double entry
        
        if(typeofs_of_betas[1] == "character"){
          x_beta_fullRE_DMSL <- rep(NA, length(list_betas[[which(typeofs_of_betas == "double")[1]]]))
        }
        
        if(typeofs_of_betas[2] == "character"){
          x_beta_diagRE_DMDL <- rep(NA, length(list_betas[[which(typeofs_of_betas == "double")[1]]]))
        }
        
        if(typeofs_of_betas[3] == "character"){
          x_beta_fullRE_M <- rep(NA, length(list_betas[[which(typeofs_of_betas == "double")[1]]]))
        }
        
        if( (length(x_beta_fullRE_DMSL) != length(x_beta_diagRE_DMDL)) | (length(x_beta_fullRE_DMSL) != length(x_beta_fullRE_M)) ){
          warning(paste0(ct, ': the number of log-ratios is not consistent'))
          return(NULL)
        }
      }
      
    }
    cbind.data.frame(rbind.data.frame(
      cbind.data.frame(fullRE_DMSL=select_slope_2(x_beta_fullRE_DMSL),
                       diagRE_DMDL=select_slope_2(x_beta_diagRE_DMDL),
                       fullRE_M=select_slope_2(x_beta_fullRE_M),
                       beta_type='Slope'),
      cbind.data.frame(fullRE_DMSL=select_intercept(x_beta_fullRE_DMSL),
                       diagRE_DMDL=select_intercept(x_beta_diagRE_DMDL),
                       fullRE_M=select_intercept(x_beta_fullRE_M),
                       beta_type='Intercept')),
      ct=ct)
  }
  ))
  
  ## if something hasn't converged, remove the value
  bad_ct_fullRE_DMSL <- c(enough_samples[sapply(enough_samples,
                                                function(ct) (typeof(model_fullRE_DMSL_list[[ct]])) == 'character')],
                          enough_samples[sapply(enough_samples,
                                                function(ct)  try(model_fullRE_DMSL_list[[ct]]$pdHess)) != "TRUE"])
  bad_ct_diagRE_DMDL <- c(enough_samples[sapply(enough_samples,
                                                function(ct) (typeof(model_diagRE_DMDL_list[[ct]])) == 'character')],
                          enough_samples[sapply(enough_samples, function(ct)  try(model_diagRE_DMDL_list[[ct]]$pdHess)) != "TRUE"])
  bad_ct_fullRE_M <- c(enough_samples[sapply(enough_samples,
                                             function(ct) (typeof(model_fullRE_M_list[[ct]])) == 'character')],
                       enough_samples[sapply(enough_samples, function(ct)  try(model_fullRE_M_list[[ct]]$pdHess)) != "TRUE"])
  
  .x$fullRE_DMSL[(.x$ct %in% bad_ct_fullRE_DMSL)] = NA
  .x$diagRE_DMDL[(.x$ct %in% bad_ct_diagRE_DMDL)] = NA
  .x$fullRE_M[(.x$ct %in% bad_ct_fullRE_M)] = NA
  
  .x$ct2=renaming_pcawg[,2][match(.x$ct, renaming_pcawg[,1])]
  
  return(.x)
  
}

L_to_cov <- function(cov_vector, d){
  warning('This function was wrong until 20220316. L is a lower triangular matrix, not a symmetrical matrix')
  L <- fill_covariance_matrix(arg_d = d, arg_entries_var = rep(1, d), arg_entries_cov = cov_vector)
  L[upper.tri(L)] <- 0
  D <- diag(L%*%t(L))
  diag((D)**(-1/2)) %*% L %*% t(L) %*% diag((D)**(-1/2))
}


comparison_randomintercepts_models <- function(model_fullRE_DMSL_list, model_diagRE_DMDL_list, model_fullRE_M_list ){
  
  .x <- do.call('rbind.data.frame', lapply(enough_samples, function(ct){
    x_RE_fullRE_DMSL <- try(python_like_select_name(model_fullRE_DMSL_list[[ct]]$par.random, "u_large"))
    x_RE_diagRE_DMDL <- try(python_like_select_name(model_diagRE_DMDL_list[[ct]]$par.random, "u_large"))
    x_RE_fullRE_M <- try(python_like_select_name(model_fullRE_M_list[[ct]]$par.random, "u_large"))
    if( (length(x_RE_fullRE_DMSL) != length(x_RE_diagRE_DMDL)) | (length(x_RE_fullRE_DMSL) != length(x_RE_fullRE_M)) ){
      ## if we don't have results for any, remove from the analysis
      list_RE <- list(x_RE_fullRE_DMSL, x_RE_diagRE_DMDL, x_RE_fullRE_M)
      typeofs_of_RE <- sapply(list_RE, typeof)
      if( all(typeofs_of_RE == "character")  ){
        return(NULL)
      }else{
        ## if we do have results for some, replace the error message by an NA string
        ## replace using the length of the first double entry
        
        if(typeofs_of_RE[1] == "character"){
          x_RE_fullRE_DMSL <- rep(NA, length(list_RE[[which(typeofs_of_RE == "double")[1]]]))
        }
        
        if(typeofs_of_RE[2] == "character"){
          x_RE_diagRE_DMSL <- rep(NA, length(list_RE[[which(typeofs_of_RE == "double")[1]]]))
        }
        
        if(typeofs_of_RE[3] == "character"){
          x_RE_fullRE_M <- rep(NA, length(list_RE[[which(typeofs_of_RE == "double")[1]]]))
        }
        
        if( (length(x_RE_fullRE_DMSL) != length(x_RE_diagRE_DMDL)) | (length(x_RE_fullRE_DMSL) != length(x_RE_fullRE_M)) ){
          warning(paste0(ct, ': the number of log-ratios is not consistent'))
          return(NULL)
        }
      }
      
    }
    
    ## put the coefficients in matrix form
    ## get the number of log-ratios, d-1
    dmin1 <- (names(table(sapply(list(model_fullRE_DMSL_list, model_diagRE_DMDL_list, model_fullRE_M_list), function(i)    as.numeric(try(length(python_like_select_name(i[[ct]]$par.fixed, 'beta'))/2))))))
    if(length(dmin1) == 1){
      ## there should only be one, shared, d-1
      dmin1 <- as.numeric(dmin1)
    }else{
      stop(paste0('Models do not agree on number of log-ratios. CT: ', ct))
    }
    
    x_RE_fullRE_DMSL <- matrix(x_RE_fullRE_DMSL, ncol=dmin1, byrow=F)
    x_RE_diagRE_DMDL <- matrix(x_RE_diagRE_DMDL, ncol=dmin1, byrow=F)
    x_RE_fullRE_M <- matrix(x_RE_fullRE_M, ncol=dmin1, byrow=F)
    
    bad_fullRE_DMSL=F
    bad_diagRE_DMDL=F
    bad_fullRE_M=F
    ## if something hasn't converged, set all the random coefficients to NA
    if((typeof(model_fullRE_DMSL_list[[ct]]) == "character") ){
      bad_fullRE_DMSL=T
    }else{
      if(try(!(model_fullRE_DMSL_list[[ct]]$pdHess))){
        bad_fullRE_DMSL=T
      }
    }
    if(bad_fullRE_DMSL){
      x_RE_fullRE_DMSL <- matrix(NA, nrow = nrow(x_RE_fullRE_DMSL), ncol=ncol(x_RE_fullRE_DMSL))
    }
    #----
    if((typeof(model_diagRE_DMDL_list[[ct]]) == "character") ){
      bad_diagRE_DMDL=T
    }else{
      if(try(!(model_diagRE_DMDL_list[[ct]]$pdHess))){
        bad_diagRE_DMDL=T
      }
    }
    if(bad_diagRE_DMDL){
      x_RE_diagRE_DMDL <- matrix(NA, nrow = nrow(x_RE_diagRE_DMDL), ncol=ncol(x_RE_diagRE_DMDL))
    }
    #-----
    if((typeof(model_fullRE_M_list[[ct]]) == "character") ){
      bad_fullRE_M=T
    }else{
      if(try(!(model_fullRE_M_list[[ct]]$pdHess))){
        bad_fullRE_M=T
      }
    }
    if(bad_fullRE_M){
      x_RE_fullRE_M <- matrix(NA, nrow = nrow(x_RE_fullRE_M), ncol=ncol(x_RE_fullRE_M))
    }
    
    ## for each patient using the x_RE_fullRE_DMSL intercepts, get the distance to the intercepts of the other two models
    dist_DMSLs <- sapply(1:nrow(x_RE_fullRE_DMSL), function(i){
      if(all(is.na(x_RE_fullRE_DMSL[i,])) | all(is.na(x_RE_diagRE_DMDL[i,]))){
        NA
      }else{
        dist(rbind(x_RE_fullRE_DMSL[i,], x_RE_diagRE_DMDL[i,]))
      }
    })
    dist_fullREs <- sapply(1:nrow(x_RE_fullRE_DMSL), function(i){
      if(all(is.na(x_RE_fullRE_DMSL[i,])) | all(is.na(x_RE_fullRE_M[i,]))){
        NA
      }else{
        dist(rbind(x_RE_fullRE_DMSL[i,], x_RE_fullRE_M[i,]))
      }
    })
    
    cbind.data.frame(melt(list(dist_DMSLs=dist_DMSLs, dist_fullREs=dist_fullREs)),
                     ct=ct)
  }
  ))
  
  .x$ct2=renaming_pcawg[,2][match(.x$ct, renaming_pcawg[,1])]
  
  return(.x)
  
}

give_first_col <- function(i){
  if(is.null(dim(i))){i[1]}else{i[,1]}
}

give_pairs_with_mvn <- function(mat_mvn, lims=NULL){
  require(car)
  if(is.null(lims)){
    dataEllipse(mat_mvn[,1], mat_mvn[,2], levels=c(0.95), xlim=c(min(mat_mvn[,1])-1,max(mat_mvn[,1])+1),
                ylim=c(min(mat_mvn[,2])-1,max(mat_mvn[,2])+1)) 
  }else{
    dataEllipse(mat_mvn[,1], mat_mvn[,2], levels=c(0.95), xlim=lims,
                ylim=lims)
  }
}

remove_all_NA <- function(i) i[!(colSums(apply(i, 1, is.na)) > 0),]

give_pairs_with_mvn_wrapper <- function(exposures_arg, zero_to_NA=F, common_lims=F){
  if(zero_to_NA){
    exposures_arg[exposures_arg == 0] <- NA
  }
  common_lims_vec <- c(min(exposures_arg, na.rm=T), max(exposures_arg, na.rm=T))
  par(mfrow=c(ncol(exposures_arg), ncol(exposures_arg)), mar=c(0,0,0,0))
  for(ii in 1:ncol(exposures_arg)){
    for(jj in 1:ncol(exposures_arg)){
      if(ii != jj){
        if(common_lims){
          give_pairs_with_mvn(remove_all_NA(exposures_arg[,c(ii, jj)]), lims=common_lims_vec)
        }else{
          give_pairs_with_mvn(remove_all_NA(exposures_arg[,c(ii, jj)]))
        }
      }else{plot.new()}
    }
  }
}

adjust_all <- function(i, new_version=F, method='BH'){
  if(!new_version){
    if(method != "holm"){
      stop('by default, in the past, using SNV, the default method has been holm')
    }
    stopifnot(colnames(i)[ncol(i)] == "true")
    cbind.data.frame(apply(i[,-ncol(i)], 2, p.adjust), true=i[,ncol(i)])
  }else{
    cbind.data.frame(apply(i[,colnames(i) != "true"], 2, p.adjust, method=method), true=i[,colnames(i) == "true"])
  }
}

