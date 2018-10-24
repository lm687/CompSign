#' Normalised signatures of 560 Breast Cancer genomes from the paper
#' Landscape of somatic mutations in 560 breast cancer whole-genome sequences
#'
#'It contains 12 signatures
#'
#' @format A merged_compositional object with a 560x12 exposure matrix and
#' a data frame with clinical data, of dimension 560x47.
#' \describe{
#'   \item{price}{price, in US dollars}
#'   \item{carat}{weight of the diamond, in carats}
#'   \item donor_gender
#'   \item donor_age_at_diagnosis
#'   \item donor_age_at_last_follow.up
#'   \item specimen_type
#'   \item donor_vital_status
#'   \item disease_status_last_follow.up
#'   \item donor_relapse_interval__in_DAYS
#'   \item donor_survival_time_in_DAYS
#'   \item donor_interval_of_last_follow.up_in_DAYS
#'   \item tumour_grade
#'   \item T_stage
#'   \item N_stage
#'   \item M_stage
#'   \item known_germline_mutations
#'   \item sample_removed_pre_or_post.treatment
#'   \item source_of_normal
#'   \item X._tumour_cellularity
#'   \item smoking.history
#'   \item other.exposure
#'   \item parity
#'   \item age_at_birth_of_first_child_.years.
#'   \item oral_contraception_exposure
#'   \item oral_contraception_years
#'   \item menopausal_status_at_diagnosis
#'   \item HRT_history_years
#'   \item number_of_positive_lymph_nodes
#'   \item DCIS_Grade
#'   \item Histopathological_subtype
#'   \item tubule_score
#'   \item pleomorphism_score
#'   \item mitotic_score
#'   \item total_mitoses
#'   \item grade
#'   \item Lymphovascular_invasion
#'   \item Lymphocyte_infiltration
#'   \item tumour_border
#'   \item Central_scar.fibrosis
#'   \item Necrosis
#'   \item X._Invasive_tumour
#'   \item X._CIS
#'   \item X._Stroma
#'   \item X._Lymphocytes
#'   \item X._Adipose_tissue
#'   \item X._Normal_epithelium
#'   \item final.ER
#'   \item final.PR
#'   \item final.HER2
#'   \item{CONTINUE}{yadayada}
#'   ...
#' }
#' @source \url{https://www.nature.com/articles/nature17676}
"Breast560"
