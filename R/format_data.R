#' Data preparation for time-varying mediaitor
#'
#' @description Internal funciton for data manipulation. Change long format to wide format.
#'
#' @param data Data set to be sued
#' @param id ID variable per subject.
#' @param trt Intervention/Exposure variable
#' @param med Name of the time-varying mediator.
#' @param y Name of the outcome variable.
#' @param cov A vector of the covariates' name.
#' @param time Time variable of the per row observed.
#'
#' @importFrom dplyr group_by arrange ungroup select mutate
#' @importFrom tidyr spread %>%
#'
#' @export
#'
#'

format_data <- function(data,
                        id,
                        trt,
                        time,
                        med,
                        y,
                        cov){

  data <- data[, c(id, trt, y, med, time, cov)]
  colnames(data)[1:5] <-c("id", "trt", "y", "med", "time")

  dt_iptw <- data %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(tm = 1:n()) %>%
    ungroup() %>%
    select(-time) %>%
    mutate(tm = paste0("med_", tm)) %>%
    spread(tm, med)
  out <- as.data.frame(dt_iptw)
  colnames(out)[1:3] <- c(id, trt, y)
  out

}
