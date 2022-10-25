## Object of class <exposures>
setClass("exposures",
         representation = representation(
           cancer_type="character",
           type_classification = "character",
           number_categories = "numeric",
           id_categories = "character",
           active_signatures = "character", ## active signatures for this cancer type
           count_matrices_all = "list", ## for each of the categories
           count_matrices_active = "list", ## for each of the categories, only active signatures
           sample_names = "character",
           is_null_active = "logical",
           is_empty = "character",
           modification = "character"
         ))



#' Function to load an exposures object
#' @param ct: cancer type, if the data are organised as in the PCAWG dataset. Otherwise, name of ile
#' @param typedata: type of data (simulation/nucleotidesubstitution1/nucleotidesubstitution3/signatures/etc.), if the data are organised as in the PCAWG dataset
#' @param simulation: boolean indicating if the data are simulations (i.e. not organised as in the PCAWG dataset)
#' @param path_to_data: path to folder which includes the file to upload
#' @param read_directly: boolean indicating if <ct> includes the filename
#' @param old_version_creating_X_Z: boolean indicating if the input file contains exposures organised as in <<Display A>>
#' @param load_all_sigs: load all signatures (otherwise, only active signatures are loaded)
#' @param override_warning_X_Z: allow X and Z to be generated assuming that the input file contains exposures organised as in <<Display A>>
#' @return the covariance matrix
load_PCAWG = function(ct, typedata, simulation=FALSE, path_to_data="../../data/", read_directly=FALSE, old_version_creating_X_Z=F, load_all_sigs=F, override_warning_X_Z=F){
  if(simulation | read_directly){
    fle = ct
    print(ct)
    cat('Reading file', fle, '\n')
  }else{
    fle = paste0(path_to_data, "roo/", ct, '_', typedata, "_ROO.RDS" )
  }
  objects_sigs_per_CT_features <- readRDS(fle)
  
  if(simulation){
    ##' get the first element of the list, because the second is the
    ##' set of the parameters used in the creation of the dataset
    objects_sigs_per_CT_features0 <- objects_sigs_per_CT_features
    objects_sigs_per_CT_features = objects_sigs_per_CT_features[[1]]
  }
  
  if(length(objects_sigs_per_CT_features) == 1){
    if(typeof(objects_sigs_per_CT_features) != "S4"){
      if(is.na(objects_sigs_per_CT_features)){
        return(NA)
      }
    }
  }
  
  if(typedata %in% c("nucleotidesubstitution1", "simulation")){
    objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
  }else if( grepl("nucleotidesubstitution3", typedata)){
    objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
    if( typedata == 'nucleotidesubstitution3MSE' ){
      objects_sigs_per_CT_features = lapply(objects_sigs_per_CT_features, t) ## need to transpose output of MSE contexts
    }
  }else if(grepl("signatures", typedata)){
    if(is.null(attr(objects_sigs_per_CT_features,"count_matrices_active")[[1]]) | (length(attr(objects_sigs_per_CT_features,"count_matrices_active")[[1]]) == 0) | load_all_sigs){
      ## no active signatures, or we want to load all signatures
      objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_all")
    }else{
      objects_sigs_per_CT_features = attr(objects_sigs_per_CT_features,"count_matrices_active")
    }
    objects_sigs_per_CT_features = lapply(objects_sigs_per_CT_features, function(i){
      rwn = rownames(i)
      .x = apply(i, 2, as.numeric)
      rownames(.x) = rwn
      round(.x)
    })
  }else{
    stop('Check <typedata> argument')
  }
  
  d = ncol(objects_sigs_per_CT_features[[1]]) ## number of signatures or features
  n = nrow(objects_sigs_per_CT_features[[1]]) ## number of samples
  
  if(typedata == "simulation"){
    ## new 20220517
    Z = objects_sigs_per_CT_features0$Z_sim
    X = objects_sigs_per_CT_features0$X_sim
  }else{
    if(!override_warning_X_Z){
      stop('Beware that this might be adding the incorrect X and Z')
    }
    if(all(rownames(objects_sigs_per_CT_features[[1]]) == rownames(objects_sigs_per_CT_features[[2]]))){
      old_version_creating_X_Z = T
    }else{
      if(typedata %in% c("signaturesMSE", "nucleotidesubstitution3MSE")){
        old_version_creating_X_Z = T
        rownames(objects_sigs_per_CT_features[[1]]) <- gsub(".consensus.20160830.somatic.snv_mnv_SUBCLONAL.vcf.gz", "", rownames(objects_sigs_per_CT_features[[1]]))
        rownames(objects_sigs_per_CT_features[[2]]) <- gsub(".consensus.20160830.somatic.snv_mnv_SUBCLONAL.vcf.gz", "", rownames(objects_sigs_per_CT_features[[2]]))
      }else{
        stop('Patients in input data are rearranged. Could not create matrices X and Z')
      }
    }
    
    if(old_version_creating_X_Z){
      ## used directly up until 2 August 2021
      
      # matrix of fixed effects
      X = matrix(NA, nrow=2, ncol=2*n)
      X[1,] = 1
      X[2,] = rep(c(0,1), each=n)
      
      # covariate matrix. two samples per patient. first all the samples from group 1, then all the samples from group 2
      Z0 = matrix(0, nrow=n, ncol=n)
      diag(Z0) = 1
      Z = t(rbind(Z0, Z0))
      
    }
  }
  ## The counts
  W = rbind(objects_sigs_per_CT_features[[1]], objects_sigs_per_CT_features[[2]]) ## used to be inside <old_version_creating_X_Z>
  
  return(list(x=t(X), z=t(Z), Y=W))
}

#' Modify the dataset, e.g. by changing BDL, transformating of the data, imputating values 
#' @param dataset_obj: dataset of type exposures_inputTMB
#' @param model_arg: name of model that will be used in TMB
#' @param imputation_arg: type of imputation (1em2, 1em4, 0, multLN, multRepl, robust)
#' @param transformation_arg: type of transformation (ALR, ILR)
#' @param BDL_arg: Below Detection Limit Threshold value
#' @param output_arg: filename of output
modify_dataset <- function(dataset_obj, model_arg, imputation_arg, transformation_arg, BDL_arg, output_arg){
  
  dataset_obj$Y[dataset_obj$Y <= BDL_arg] = 0
  
  if(grepl('bernoulli', model_arg)){
    cat('Not modifying imput. We are analysing presence/absence of signatures\n')
    
  }else{
    
    ## analysing the exposures
    if(sum(dataset_obj$Y == 0) == 0){
      cat('No zero instances\n')
    }
    
    if(imputation_arg == '1em2'){
      if(sum(dataset_obj$Y == 0) == 0){
        cat('No zero intances\n')
        saveRDS(object = NA, file = output_arg)
        quit()
      }
      dataset_obj$Y <- impute(dataset_obj$Y, 1e-2)
    }else if(imputation_arg == '1em4'){
      if(sum(dataset_obj$Y == 0) == 0){
        cat('No zero intances\n')
        saveRDS(object = NA, file = output_arg)
        quit()
      }
      dataset_obj$Y <- impute(dataset_obj$Y, 1e-4)
    }else if(imputation_arg == '0'){
      print('No imputation')
    }else if(imputation_arg == 'multLN'){
      if(sum(dataset_obj$Y == 0) == 0){
        cat('No zero intances\n')
        saveRDS(object = NA, file = output_arg)
        quit()
      }
      dataset_obj$Y <- zCompositions::multLN(dataset_obj$Y, label = 0, dl = rep(BDL_arg, ncol(dataset_obj$Y)))
    }else if(imputation_arg == 'multRepl'){
      if(sum(dataset_obj$Y == 0) == 0){
        cat('No zero intances\n')
        saveRDS(object = NA, file = output_arg)
        quit()
      }
      dataset_obj$Y <- zCompositions::multRepl(dataset_obj$Y, label = 0, dl = rep(BDL_arg, ncol(dataset_obj$Y)))
    }else if(imputation_arg == 'robust'){
      if(transformation_arg == 'ALR'){
        dataset_obj$Y <- as(compositions::alr(dataset_obj$Y), 'matrix')
      }else if(transformation_arg == 'ILR'){
        dataset_obj$Y <- as(compositions::ilr(dataset_obj$Y), 'matrix')
      }else{
        stop('Robust imputation is only allowed for ALR and ILR transformations')
      }
    }else{
      stop('Specify appropriate imputation')
    }
    
    ## if we are using the robust imputation we already have the transformed values above
    if(imputation_arg != 'robust'){
      if(transformation_arg == 'ALR'){
        dataset_obj$Y <- as(compositions::alr(dataset_obj$Y), 'matrix')
      }else if(transformation_arg == 'ILR'){
        if(grepl('partialILR', model_arg)){ ## changed 20220223. it used to be the case that for partialILR without imputation the robust ILR would be used
          cat('Removing runs with zeros\n')
          dataset_obj$Y <- give_partial_irl(dataset_obj$Y)
          .keep = !(rowSums(dataset_obj$Y) == 0)
          dataset_obj$Y = dataset_obj$Y[.keep,]
          dataset_obj$x = dataset_obj$x[.keep,]
        }else{
          dataset_obj$Y <- as(compositions::ilr(dataset_obj$Y), 'matrix')
        }
      }else{
        stop('Specify appropriate transformation')
      }
    }
  }
  
  ## modify matrix z if necessary, if any samples have been removed
  if(exists(".keep")){
    dataset_obj$z = dataset_obj$z[.keep,]
    dataset_obj$z <- dataset_obj$z[,colSums(dataset_obj$z) > 0]
  }
  
  return(dataset_obj)
  
}