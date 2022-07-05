
#library(boot)

## Bootstrap helpers

print_CI <- function(boot.ci_output, stat_name=NULL){
  ## Prints boot::boot.ci function call output in a formated manner
  cat("Confidence Intervals based on",boot.ci_output$R,"bootstrap replicates")
  if(!is.null(stat_name)){
    cat(" for", stat_name, sep=" ")
  }
  cat(":\n")
  if(!is.null(boot.ci_output[["normal"]])){
    cat("\t Normal  (", boot.ci_output$normal[[1]]*100 ,"%):\t(", boot.ci_output$normal[[2]], ",", boot.ci_output$normal[[3]], ")\n")
  }
  if(!is.null(boot.ci_output[["basic"]])){
    cat("\t Basic  (", boot.ci_output$basic[[1]]*100 ,"%):\t(", boot.ci_output$basic[[4]], ",", boot.ci_output$basic[[5]], ")\n")
  }
  if(!is.null(boot.ci_output[["student"]])){
    cat("\t Student (", boot.ci_output$student[[1]]*100 ,"%):\t(", boot.ci_output$student[[4]], ",", boot.ci_output$student[[5]], ")\n")
  }
  if(!is.null(boot.ci_output[["percent"]])){
    cat("\t Percentile (", boot.ci_output$percent[[1]]*100 ,"%):\t(", boot.ci_output$percent[[4]], ",", boot.ci_output$percent[[5]], ")\n")
  }
  if(!is.null(boot.ci_output[["bca"]])){
    cat("\t BCA  (", boot.ci_output$bca[[1]]*100 ,"%):\t(", boot.ci_output$bca[[4]], ",", boot.ci_output$bca[[5]], ")\n")
  }
}


## Functions to obtain the statistics from the data, to be used with boot package
CTC_boot_stats <- function(data, indices, k=floor(nrow(data)^0.4), both_tails=FALSE) {
  d <- data[indices,] # allows boot to select sample
  coeff12 <- causal_tail_coeff(d[[1]],d[[2]], k=k, both_tails=both_tails)
  coeff21 <- causal_tail_coeff(d[[2]],d[[1]], k=k, both_tails=both_tails)
  return(c(coeff12 = coeff12, coeff21 = coeff21, diff12 = coeff12 - coeff21))
}
LGPDCTC_boot_stats <- function(data, indices, k=floor(nrow(data)^0.4), constrained_fit=FALSE, threshold_q = 0.95, parametric_F1=TRUE, both_tails = FALSE) {
  d <- data[indices,] # allows boot to select sample
  if(constrained_fit){
    coeff12 <- constrained_LGPD_CTC(d[[1]], d[[2]], d[[3]], k=k, threshold_q=threshold_q, parametric_F1=parametric_F1, both_tails=both_tails)
    coeff21 <- constrained_LGPD_CTC(d[[2]], d[[1]], d[[3]], k=k, threshold_q=threshold_q, parametric_F1=parametric_F1, both_tails=both_tails)
  }else{
    coeff12 <- LGPD_causal_tail_coeff(d[[1]], d[[2]], d[[3]], k=k, threshold_q=threshold_q, parametric_F1=parametric_F1, both_tails=both_tails)
    coeff21 <- LGPD_causal_tail_coeff(d[[2]], d[[1]], d[[3]], k=k, threshold_q=threshold_q, parametric_F1=parametric_F1, both_tails=both_tails)
  }
  return(c(coeff12 = coeff12, coeff21 = coeff21, diff12 = coeff12 - coeff21))
}
expGPDCTC_boot_stats <- function(data, indices, k=floor(nrow(data)^0.4), threshold_q = 0.95, parametric_F1=TRUE, both_tails = FALSE, constrained_fit=NULL) {
  d <- data[indices,] # allows boot to select sample
  coeff12 <- expGPD_causal_tail_coeff(d[[1]], d[[2]], d[[3]], k=k, threshold_q=threshold_q, parametric_F1=parametric_F1, both_tails=both_tails)
  coeff21 <- expGPD_causal_tail_coeff(d[[2]], d[[1]], d[[3]], k=k, threshold_q=threshold_q, parametric_F1=parametric_F1, both_tails=both_tails)
  return(c(coeff12 = coeff12, coeff21 = coeff21, diff12 = coeff12 - coeff21))
}

both_coeffs <- function(data, indices, k=floor(nrow(data)^0.4), both_tails=TRUE) {
  #Deprecated
  d <- data[indices,] # allows boot to select sample
  coeff12 <- causal_tail_coeff(d[[1]],d[[2]], k=k, both_tails=both_tails)
  coeff21 <- causal_tail_coeff(d[[2]],d[[1]], k=k, both_tails=both_tails)
  return(c(coeff12 = coeff12, coeff21 = coeff21))
}

diff_coeffs <- function(data, indices, k=floor(nrow(data)^0.4), both_tails=TRUE) {
  #Deprecated
  d <- data[indices,] # allows boot to select sample
  coeff12 <- causal_tail_coeff(d[[1]],d[[2]], k=k, both_tails=both_tails)
  coeff21 <- causal_tail_coeff(d[[2]],d[[1]], k=k, both_tails=both_tails)
  return(coeff12 - coeff21)
}
