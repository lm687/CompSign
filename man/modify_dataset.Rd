% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROO_functions.R
\name{modify_dataset}
\alias{modify_dataset}
\title{Modify the dataset, e.g. by changing BDL, transformating of the data, imputating values}
\usage{
modify_dataset(
  dataset_obj,
  model_arg,
  imputation_arg,
  transformation_arg,
  BDL_arg,
  output_arg
)
}
\arguments{
\item{dataset_obj:}{dataset of type exposures_inputTMB}

\item{model_arg:}{name of model that will be used in TMB}

\item{imputation_arg:}{type of imputation (1em2, 1em4, 0, multLN, multRepl, robust)}

\item{transformation_arg:}{type of transformation (ALR, ILR)}

\item{BDL_arg:}{Below Detection Limit Threshold value}

\item{output_arg:}{filename of output}
}
\description{
Modify the dataset, e.g. by changing BDL, transformating of the data, imputating values
}
