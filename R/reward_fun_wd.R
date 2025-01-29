reward_fun_wd <- function(episimdata,alpha,alpha_d,ovp,dovp,C_target,C_target_pen,D_target,D_target_pen,actions,ii,jj) {

  C_err_pred <- abs(episimdata[ii, 'C'] - C_target)
  R_err_pred <- abs(episimdata[ii, 'Re'] - R_target)
  D_err_pred <- abs(episimdata[ii, 'Deaths'] - D_target)
  over_pen <- 0.0
  over_pend_d <- 0.0

  if (episimdata[ii, 'C'] > C_target_pen) {
    over_pen <- ovp
  }

  if (episimdata[ii, 'Deaths'] > D_target_pen) {
    over_pen_d <- dovp
  }

  rew <- -alpha * C_err_pred - over_pen - actions[jj,"cost_of_NPI"] -alpha_d * D_err_pred - over_pend_d

  return(rew)
}
