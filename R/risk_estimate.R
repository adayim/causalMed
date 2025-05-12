

# Risk difference and risk ratio calculation function
risk_estimate1 <- function(data_list, ref_int, intervention, return_data) {
  ref_dat <- data_list[[ref_int]]
  vs_nam <- setdiff(names(intervention), ref_int)
  out <- sapply(vs_nam, function(x) {
    tmp_dt <- data_list[[x]]
    if(return_data){
      ref_mean <- sum(unlist(ref_dat[["Pred_Y"]])) / length(ref_dat[["Pred_Y"]])
      vs_mean <- sum(unlist(tmp_dt[["Pred_Y"]])) / length(tmp_dt[["Pred_Y"]])
    }else{
      ref_mean <- ref_dat
      vs_mean <- tmp_dt
    }
    data.table(
      Intervention = c(paste(x, ref_int, sep = " - "),
                       paste(x, ref_int, sep = " / ")),
      Risk_type = c("Difference", "Ratio"),
      Estimate = c(vs_mean - ref_mean, vs_mean / ref_mean)
    )
  }, simplify = FALSE)
  data.table::rbindlist(out, use.names = TRUE)
}

risk_estimate2 <- function(data_list, ref_int, intervention, return_data) {
  ref_dat <- data_list[data_list$Intervention==ref_int,]
  vs_nam <- setdiff(names(intervention), ref_int)
  out <- sapply(vs_nam, function(x) {
    tmp_dt <- data_list[data_list$Intervention==x,]
    if(return_data){
      ref_mean <- sum(ref_dat[["Est"]]) / length(ref_dat[["Est"]])
      vs_mean <- sum(tmp_dt[["Est"]]) / length(tmp_dt[["Est"]])
    }else{
      ref_mean <- ref_dat$Est
      vs_mean <- tmp_dt$Est
    }
    data.table(
      Intervention = c(paste(x, ref_int, sep = " - "),
                       paste(x, ref_int, sep = " / ")),
      Risk_type = c("Difference", "Ratio"),
      Estimate = c(vs_mean - ref_mean, vs_mean / ref_mean)
    )
  }, simplify = FALSE)
  data.table::rbindlist(out, use.names = TRUE)
}

# Calculate mediation effect for point estimates
risk_estimate.mediation1 <- function(data_list, return_data) {
  if(return_data){
    phi_11 <- sum(data_list$Ph11[["Pred_Y"]]) / length(data_list$Ph11[["Pred_Y"]])
    phi_00 <- sum(data_list$Phi00[["Pred_Y"]]) / length(data_list$Phi00[["Pred_Y"]])
    phi_10 <- sum(data_list$Phi10[["Pred_Y"]]) / length(data_list$Phi10[["Pred_Y"]])
    data.table(Effect = c("Indirect effect", "Direct effect", "Total effect"),
               Est = c(phi_11 - phi_10, phi_10 - phi_00, phi_11 - phi_00)
    )
  }else{
    data.table(Effect = c("Indirect effect", "Direct effect", "Total effect"),
               Est = c(data_list$Ph11 - data_list$Phi10,
                       data_list$Phi10 - data_list$Phi00,
                       data_list$Ph11 - data_list$Phi00)
    )
  }
}

# Calculate mediation effect for bootstrap
risk_estimate.mediation2 <- function(data_list, return_data) {
  if(return_data){
    phi_11 <- sum(data_list$Ph11) / length(data_list$Ph11)
    phi_00 <- sum(data_list$Phi00) / length(data_list$Phi00)
    phi_10 <- sum(data_list$Phi10) / length(data_list$Phi10)
    data.table(Effect = c("Indirect effect", "Direct effect", "Total effect"),
               Est = c(phi_11 - phi_10, phi_10 - phi_00, phi_11 - phi_00)
    )
  }else{
    data.table(Effect = c("Indirect effect", "Direct effect", "Total effect"),
               Est = c(data_list$Ph11 - data_list$Phi10,
                       data_list$Phi10 - data_list$Phi00,
                       data_list$Ph11 - data_list$Phi00)
    )
  }
}

