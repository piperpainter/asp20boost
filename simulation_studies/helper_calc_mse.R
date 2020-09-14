calc_MSE <- function(true_beta, true_gamma, est_beta, est_gamma) {
  result_list <- list()
  result_list$beta <- 0
  result_list$gamma <- 0
  result_list$pooled <- 0

  for(i in 1:2){
    result_list$beta <- result_list$beta + (true_beta[i] - est_beta[i])^2
    result_list$gamma <- result_list$gamma + (true_gamma[i] - est_gamma[i])^2
  }
  result_list$pooled <- result_list$beta + result_list$gamma
  return(result_list)
}
