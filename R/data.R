#' @format ## `PancEndocrine_signaturesMSE`
#' A list of length three:
#' \describe{
#'   \item{Y}{Exposures of 140 samples (70 patients) and 13 signatures}
#'   \item{x}{Match of samples and groups (second column), and including the intercept (first column)}
#'   \item{z}{Match of patients and samples}
#' }
#' @source https://dcc.icgc.org/pcawg
"PancEndocrine_signaturesMSE"


#' @format ## `simplified_object`
#' A list of length three:
#' \describe{
#'   \item{Y}{Exposures of 34 samples (17 patients) and 6 nucleotide changes}
#'   \item{x}{Match of samples and groups (second column), and including the intercept (first column), as a binary matrix}
#'   \item{z}{Match of patients and samples, as a binary matrix}
#' }
#' @source https://dcc.icgc.org/pcawg
"simplified_object"

#' @format ## `ProstAdenoCA_chrom`
#' A list of length two:
#' \describe{
#'   \item{all_exposures_ct}{List of length 208 (the number of patients in the PCAWG ProstAdenoCA cohort) matrices, where each matrix is a signature x chromosome matrix of signature exposures, for all signatures.}
#'   \item{all_exposures_ct_active}{List of length 208 (the number of patients in the PCAWG ProstAdenoCA cohort) matrices, where each matrix is a signature x chromosome matrix of signature exposures, for active signatures in this cancer type.}
#' }
#' @source https://dcc.icgc.org/pcawg
"ProstAdenoCA_chrom"
        
#' @format ## `obj_multilambda`
#' Simulated data for running a model with patient-specific overdispersion parameters, as shown in the vignette. This requires multiple observations per patient. A list of length three:
#' \describe{
#'   \item{Y}{Exposures of 1200 samples (40 patients) and 5 mutational categories}
#'   \item{x}{Match of samples and groups (second column), and including the intercept (first column), as a binary matrix}
#'   \item{z}{Match of patients and samples, as a binary matrix}
#' }
"obj_multilambda"

#' @format ## `obj_multilambda_parameters`
#' Parameters used to simulate data in obj_multilambda, as shown in the vignette.
#' }
"obj_multilambda_parameters"

#' @format ## `res_patient_lambda`
#'Results from the multilambda model, as shown in the vignette.
#' }
"res_patient_lambda"


#' @format ## `ProstAdenoCA_trinucleotide`
#' Object that can be used as input in signature extraction functions (see e.g. extract_sigs_TMB_obj). A list of length three:
#' \describe{
#'   \item{Y}{Exposures of 416 samples (208 patients) and 96 trinucleotide changes, from the PCAWG ProstAdenoCA cohort.}
#'   \item{x}{Match of samples and groups (second column), and including the intercept (first column), as a binary matrix}
#'   \item{z}{Match of patients and samples, as a binary matrix}
#' }
#' @source https://dcc.icgc.org/pcawg
"ProstAdenoCA_trinucleotide"
