
#' Simulate Data
#'
#'
#'
#' @param data Data to be used for the data generation
#' @param models Model list passed from \code{\link{Gformula}}.
#' @param intervention A vector, intervention treatment per time.
#' @param mediation_type Type of the mediation anlaysis.
#'
simulate_data <- function(data,
                          models,
                          exposure,
                          intervention = NULL,
                          mediation_type = c(NULL, "N", "I")){

  mediation_type <- match.arg(mediation_type)

  # Get the variable name of exposure
  exp_flag <- sapply(models, function(mods) as.numeric(mods$mod_type == "exposure"))
  exp_flag <- which(exp_flag == 1)
  exposure <- all.vars(formula(models[[exp_flag]]$fitted)[[2]])

  # If Intervention is given, set the treatment to given value
  if(!is.null(intervention) & is.null(mediation_type)){
    data[[exposure]] <- intervention
  }

  # if the mediation type is defined, then set the intervention to 1. But 0
  # for mediator. This is to calculate the phi_10
  if(!is.null(mediation_type)){
    data[[exposure]] <- 1
    interv0 <- parse(text = paste0(exposure, " = 0"))
  }

  # Loop through models
  for(indx in seq_along(models)){

    # Get the name of the response variable for the current model
    resp_var <- models$rsp_vars
    mod_type <- models$mod_type

    # Skip if the response variable is treatment and the intervention is defined.
    if(resp_var == exposure & !is.null(intervention))
      next

    # Perform the recode in the model
    if(!is.null(models[[indx]]$recode)){
      for(i in seq_along(models[[indx]]$recode)){
        data <- within(data, eval(parse(text = models[[indx]]$recode[i])))
      }
    }

    # If condition has been defined, apply it
    if(!is.null(models[[indx]]$subset)){
      cond <- with(data, eval(models[[indx]]$subset))
      cond <- cond & !is.na(cond) # Avoid NA in the condition list.
    }else{
      cond <- rep(TRUE, nrow(data))
    }

    if(sum(cond) != 0){
      if(mod_type == "mediator" & !is.null(mediation_type)){
        # Set the mediator's intervention to 0
        data[[resp_var]][cond] <- sim_value(within(data[cond,  ],
                                                    eval(interv0)),
                                             models[[indx]])

      }else{
        data[[resp_var]][cond] <- sim_value(data[cond, ], models[[indx]])
      }
    }

  }

  return(data)
}
