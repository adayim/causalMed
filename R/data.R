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


#' Yamamuro et al. (2021) Multi-Mediator Simulation Dataset
#'
#' A single dataset simulated from the data-generating process of the
#' Yamamuro et al. (2021) multi-mediator simulation study: a time-varying
#' binary treatment, a time-varying confounder, two sequential continuous
#' mediators, and a discrete-time survival outcome over three visits. The
#' large-sample true interventional effects for this process are documented
#' below, so estimates obtained with \code{\link{mediation}} can be compared
#' against known values.
#'
#' @format A data frame with 29482 rows and 15 variables (long at-risk
#'   format: one row per subject per visit until the event):
#' \describe{
#'   \item{id}{Unique subject identifier.}
#'   \item{time}{Visit index (0, 1, 2).}
#'   \item{V}{Time-fixed binary baseline covariate.}
#'   \item{A}{Time-varying binary treatment.}
#'   \item{L}{Time-varying continuous confounder.}
#'   \item{M1}{First mediator (continuous), affected by \code{A} and \code{L}.}
#'   \item{M2}{Second mediator (continuous), affected by \code{A}, \code{L},
#'     and \code{M1}.}
#'   \item{Y}{Event indicator (1 = event at this visit).}
#'   \item{lag1_A, lag1_L, lag1_M1, lag1_M2}{Previous-visit values (0 at
#'     baseline).}
#'   \item{L0base, M10base, M20base}{Baseline (visit 0) values of \code{L},
#'     \code{M1}, \code{M2}, carried as time-fixed columns so the Monte Carlo
#'     simulation can be seeded via \code{init_recode}.}
#' }
#'
#' @details
#' The within-visit ordering is \strong{A \eqn{\to} L \eqn{\to} M1 \eqn{\to}
#' M2 \eqn{\to} Y}; the full data-generating equations are translated from the
#' SAS \code{\%simdata} macro in the paper's supplementary material and are
#' reproduced in \code{data-raw/yamamurodata.R} in the package source
#' repository.
#'
#' The published study design generates 1000 replicate datasets of
#' \eqn{n = 1000} subjects and averages the estimates. A single replicate of
#' that size is dominated by sampling error (between-replicate SD of the
#' total effect is about 1.7 percentage points), so the dataset shipped here
#' is one draw of \eqn{n = 10000} subjects (seed 2468) from the identical
#' process, making single-dataset estimates informative.
#'
#' \strong{Large-sample true values} (risk differences in percentage points,
#' computed from the data-generating process at \eqn{n = 10^7}; Monte Carlo
#' SE in parentheses):
#' \tabular{lr}{
#'   Total effect (TE = \eqn{E[Y_1] - E[Y_0]}) \tab \eqn{-6.36} (0.011) \cr
#'   Interventional direct effect (IDE) \tab \eqn{-3.20} (0.013) \cr
#'   Interventional indirect effect via M1 \tab \eqn{-2.29} (0.011) \cr
#'   Interventional indirect effect via M2 \tab \eqn{-0.97} (0.009) \cr
#'   Mediated-interaction residual (TE \eqn{-} overall) \tab \eqn{0.10} (0.016) \cr
#' }
#' Every interventional arm, including the references, draws each mediator
#' jointly over time from an independently permuted marginal pool
#' (Yamamuro et al. 2021, Eq. 2 and Fig. 3); the total effect is the
#' natural-course contrast, so the residual row is non-zero.
#'
#' @source Simulated from the data-generating process of Yamamuro, S.,
#'   Shinozaki, T., Iimuro, S., & Matsuyama, Y. (2021). Mediational g-formula
#'   for time-varying treatment and repeated-measured multiple mediators.
#'   \emph{Statistical Methods in Medical Research}, 30(8), 1782-1799.
#'   \doi{10.1177/09622802211025988}
"yamamurodata"
