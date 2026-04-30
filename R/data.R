#' Lipid test results of 500 persons.
#'
#' A dataset containing the lipid profile and other attributes.
#'
#' @format A data frame with 2114 rows and 13 variables:
#' \describe{
#'   \item{id}{unique id of subject}
#'   \item{age0}{baseline age}
#'   \item{gender}{gender of subject}
#'   \item{class}{class membership of per-subject}
#'   \item{smoke}{Ever smoke, 0=No, 1=Yes}
#'   \item{bmi}{BMI, in kg/m^2}
#'   \item{hdl}{high density lipoprotein, mmol/L, time-varying}
#'   \item{ldl}{low density lipoprotein, mmol/L, time-varying}
#'   \item{tg}{triglyceride, mmol/L, time-varying}
#'   \item{time}{Followup time since baseline (time = 0)}
#'   \item{cvd}{cardiovascular diseases event}
#'   \item{os}{Overall Survival of cardiovascular diseases}
#' }
#' @source \url{http://www.diamondse.info/}
"lipdat"


#' Bone marrow transplant data
#'
#' The study population arose from a multicenter trial of leukemia patients and
#' comprises 137 individuals prepared for bone marrow transplants under a
#' radiation-free regimen at four medical centers.
#'
#' @format A data frame with 108714 rows and 29 variables:
#' \describe{
#'   \item{id}{unique id of subject}
#'   \item{age}{baseline age}
#'   \item{male}{male indicator}
#'   \item{cmv}{cytomegalovirus immune status (yes or no)}
#'   \item{all}{leukemia type (acute lymphocytic = 1 or acute myeloid leukemia = 0)}
#'   \item{wait}{wait time from leukemia diagnosis to transplantation}
#'   \item{day}{days since bone marrow transplant}
#'   \item{d}{indicator of death (1= yes, 0=no)}
#'   \item{gvhd}{indicator of GvHD (1= yes, 0=no)}
#'   \item{daysgvhd}{number of days since onset of GvHD}
#'   \item{daysnogvhd}{number of days without GvHD}
#'   \item{gvhdm1}{lagged GvHD}
#'   \item{relapse}{relapse indicator}
#'   \item{relapsem1}{lagged relapse indicator}
#'   \item{daysrelapse}{number of days since relapse}
#'   \item{daysnorelapse}{number of days relapse-free}
#'   \item{platnorm}{normal platelet levels (1=patient has relapsed or reached
#'   normal platelet count, 0=not in relapse or below normal platelets)}
#'   \item{platnormm1}{lagged normal platelet levels}
#'   \item{daysplatnorm}{time spent reaching normal platelet levels}
#'   \item{daysnoplatnorm}{time spent without reaching normal platelet levels}
#'   \item{censlost}{indicator of censoring due to loss-to-follow-up (1=yes, 0=no)}
#' }
#' @references Keil et al. (2014) Epidemiology, 25(6), 889-897
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4310506/}{PubMed})
"gvhd"


#' Example Dataset for a Non-Survival Outcome
#'
#' A simulated dataset with time-varying and baseline variables for 1000 subjects
#' over 5 time points, including exposure, mediator, confounders, and outcome.
#'
#' @format A data frame with 5000 rows and 13 variables:
#' \describe{
#'   \item{id}{Unique subject identifier.}
#'   \item{time}{Time variable (0 to 4).}
#'   \item{V}{Time-fixed baseline covariate.}
#'   \item{L1}{Time-varying confounder 1 (continuous).}
#'   \item{L2}{Time-varying confounder 2 (binary).}
#'   \item{A}{Time-varying binary exposure.}
#'   \item{M}{Time-varying mediator.}
#'   \item{Y_bin}{Binary outcome observed at each time point.}
#'   \item{Y_cont}{Continuous outcome observed at each time point.}
#'   \item{lag1_A}{Lagged exposure (A at previous time point).}
#'   \item{lag1_L1}{Lagged confounder L1.}
#'   \item{lag1_L2}{Lagged confounder L2.}
#'   \item{lag1_M}{Lagged mediator.}
#' }
#' 
#' @details
#' The simulated longitudinal data-generating structure can be summarized as:
#' \deqn{
#' A_t \leftarrow V, L1_{t-1}, L2_{t-1}, A_{t-1}, t;\quad
#' L1_t \leftarrow V, A_t, L1_{t-1}, t;\quad
#' L2_t \leftarrow V, A_t, L2_{t-1}, t;\quad
#' M_t \leftarrow V, A_t, L1_t, L2_t, M_{t-1}, t;\quad
#' Y_t \leftarrow V, A_t, M_t, L1_t, L2_t, A_t*M_t.
#' }
#' The same outcome model structure is used for both \code{Y_bin} and
#' \code{Y_cont}, with the appropriate outcome distribution specified for each
#' outcome type.
#'
#' @source Simulated data generated for package examples.
"nonsurvivaldata"


#' Example Dataset for a Survival Outcome
#'
#' A simulated dataset with time-varying and baseline variables for subjects
#' with a survival outcome, suitable for use with the g-formula and mediation
#' functions under survival settings.
#'
#' @format A data frame with 7113 rows and 10 variables:
#' \describe{
#'   \item{id}{Unique subject identifier.}
#'   \item{time}{Time variable (integer, starting at 0).}
#'   \item{V}{Time-fixed baseline covariate.}
#'   \item{L}{Time-varying confounder.}
#'   \item{A}{Time-varying binary exposure.}
#'   \item{M}{Time-varying mediator.}
#'   \item{Y}{Survival outcome indicator (1 = event, 0 = alive/censored).}
#'   \item{lag1_A}{Lagged exposure (A at previous time point).}
#'   \item{lag1_M}{Lagged mediator.}
#'   \item{lag1_L}{Lagged confounder.}
#' }
#' 
#' @details
#' The data-generating structure can be summarized as:
#' \deqn{
#' A_t \leftarrow V, L_{t-1}, A_{t-1}, t;\quad
#' L_t \leftarrow V, A_t, L_{t-1}, t;\quad
#' M_t \leftarrow V, A_t, L_t, M_{t-1}, t;\quad
#' Y_t \leftarrow V, A_t, M_t, L_t, A_t*M_t, t.
#' }
#'
#' Events and follow-up observations are generated only while subjects remain at
#' risk. Therefore, once \code{Y} becomes 1 at time \eqn{t}, no observations are
#' retained for that subject at subsequent time points \eqn{t+1, t+2, \ldots}.
#' 
#' @source Simulated data generated for package examples.
"survivaldata"
