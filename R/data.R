#' Bone marrow transplant data (graft-versus-host disease)
#'
#' Person-day data for the 137 leukemia patients who received bone marrow
#' transplants under a radiation-free regimen at four medical centers, as used in
#' the parametric g-formula illustration of Keil et al. (2014). Each subject
#' contributes one row per day of follow-up (day 1 up to death, loss to
#' follow-up, or day 1825), giving 108,714 rows. The dataset reproduces the
#' person-day file created in Appendix 1 of that paper: 80 deaths and 39 losses
#' to follow-up.
#'
#' GvHD (\code{gvhd}) is the exposure, death (\code{d}) is the outcome, and
#' relapse (\code{relapse}) and platelet recovery (\code{platnorm}) are
#' time-varying confounders. The three time-varying states are \emph{absorbing}
#' (once 1, they remain 1); their \code{daysX}/\code{daysnoX} counters and
#' one-day lags (\code{Xm1}) summarize the history entering each model.
#'
#' @format A data frame with 108,714 rows and 21 variables:
#' \describe{
#'   \item{id}{unique subject identifier (1-137)}
#'   \item{age}{baseline age in years}
#'   \item{male}{male indicator (1 = male, 0 = female)}
#'   \item{cmv}{baseline cytomegalovirus immune status (1 = positive, 0 = negative)}
#'   \item{all}{leukemia type (1 = acute lymphocytic, 0 = acute myeloid)}
#'   \item{wait}{waiting time from leukemia diagnosis to transplantation, in
#'   months (waiting days / 30.5)}
#'   \item{day}{days since bone marrow transplant (1-1825)}
#'   \item{d}{indicator of death at the end of day \code{day} (1 = yes, 0 = no)}
#'   \item{gvhd}{indicator of graft-versus-host disease having occurred by
#'   day \code{day} (1 = yes, 0 = no); the exposure}
#'   \item{daysgvhd}{cumulative number of days with GvHD}
#'   \item{daysnogvhd}{cumulative number of days without GvHD}
#'   \item{gvhdm1}{one-day lag of \code{gvhd}}
#'   \item{relapse}{indicator of leukemia relapse having occurred by day
#'   \code{day} (1 = yes, 0 = no)}
#'   \item{relapsem1}{one-day lag of \code{relapse}}
#'   \item{daysrelapse}{cumulative number of days since relapse}
#'   \item{daysnorelapse}{cumulative number of days relapse-free}
#'   \item{platnorm}{indicator of platelet recovery having occurred by day
#'   \code{day} (1 = platelets returned to normal levels, 0 = not yet)}
#'   \item{platnormm1}{one-day lag of \code{platnorm}}
#'   \item{daysplatnorm}{cumulative number of days with normal platelet levels}
#'   \item{daysnoplatnorm}{cumulative number of days without normal platelet
#'   levels}
#'   \item{censlost}{indicator of censoring due to loss to follow-up at day
#'   \code{day} (1 = yes, 0 = no)}
#' }
#'
#' @details
#' A worked total-effect analysis on this dataset (counterfactual mortality risk
#' had GvHD never occurred versus the natural course), mirroring the model
#' specifications in Appendix 2 of Keil et al. (2014), is walked through in
#' \code{vignette("causalMed-03-gformula")}.
#'
#' @references Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., &
#'   Cole, S. R. (2014). The parametric g-formula for time-to-event data:
#'   intuition and a worked example. \emph{Epidemiology}, 25(6), 889-897.
#'   \doi{10.1097/EDE.0000000000000160}
"gvhd"


#' Example Dataset for a Non-Survival Outcome
#'
#' A simulated dataset with time-varying and baseline variables for 3000
#' subjects over 5 time points, including exposure, mediator, confounders, and
#' an end-of-follow-up outcome. Every subject contributes all 5 rows (there is
#' no attrition), which makes this the simpler of the two simulated examples.
#'
#' @format A data frame with 15,000 rows (3000 subjects x 5 time points) and 13
#'   variables:
#' \describe{
#'   \item{id}{Unique subject identifier (1-3000).}
#'   \item{time}{Time variable (0 to 4).}
#'   \item{V}{Time-fixed baseline covariate.}
#'   \item{L1}{Time-varying confounder 1 (continuous).}
#'   \item{L2}{Time-varying confounder 2 (binary).}
#'   \item{A}{Time-varying binary exposure.}
#'   \item{M}{Time-varying mediator (continuous).}
#'   \item{Y_bin}{Binary end-of-follow-up outcome. Recorded \strong{only at
#'     \code{time == 4}}; \code{NA} at \code{time} 0-3.}
#'   \item{Y_cont}{Continuous end-of-follow-up outcome, on the same schedule as
#'     \code{Y_bin}: recorded \strong{only at \code{time == 4}}, \code{NA}
#'     otherwise.}
#'   \item{lag1_A}{Exposure at the previous time point; \code{NA} at
#'     \code{time == 0}.}
#'   \item{lag1_L1}{Previous value of \code{L1}; \code{NA} at \code{time == 0}.}
#'   \item{lag1_L2}{Previous value of \code{L2}; \code{NA} at \code{time == 0}.}
#'   \item{lag1_M}{Previous value of \code{M}; \code{NA} at \code{time == 0}.}
#' }
#'
#' @details
#' The simulated longitudinal data-generating structure can be summarized as:
#' \deqn{
#' A_t \leftarrow V, L1_{t-1}, L2_{t-1}, A_{t-1}, t;\quad
#' L1_t \leftarrow V, A_t, L1_{t-1}, t;\quad
#' L2_t \leftarrow V, A_t, L2_{t-1}, t;\quad
#' M_t \leftarrow V, A_t, L1_t, L2_t, M_{t-1}, t;\quad
#' Y \leftarrow V, A_4, M_4, L1_4, L2_4, A_4*M_4.
#' }
#' The exposure, mediator and confounders evolve at every time point, but the
#' outcome is realised once, at the end of follow-up. The same outcome
#' structure is used for both \code{Y_bin} and \code{Y_cont}, with the
#' appropriate distribution specified for each.
#'
#' Two consequences for model specification:
#' \itemize{
#'   \item Because the outcome exists at a single time point, \code{time} is
#'     constant among the rows used to fit the outcome model and must
#'     \strong{not} appear in its formula: the term is not estimable and would
#'     be dropped (with a warning) from the simulation.
#'   \item The \code{lag1_*} columns are \code{NA} at \code{time == 0}, so a
#'     baseline value must be supplied through \code{init_recode}, e.g.
#'     \code{init_recode = recodes(lag1_A = 0, lag1_L1 = 0)}. (The lag columns
#'     of \code{\link{survivaldata}} are 0-filled at baseline instead, so no
#'     initialisation is strictly required there.)
#' }
#'
#' @seealso \code{\link{survivaldata}} for the survival-outcome counterpart.
#' @source Simulated data generated for package examples.
"nonsurvivaldata"


#' Example Dataset for a Survival Outcome
#'
#' A simulated dataset with time-varying and baseline variables for 3000
#' subjects with a discrete-time survival outcome, suitable for use with the
#' g-formula and mediation functions under survival settings.
#'
#' @format A data frame with 7113 rows and 10 variables, in long at-risk format
#'   (one row per subject per time point while still at risk; 3000 subjects,
#'   2799 events, no missing values):
#' \describe{
#'   \item{id}{Unique subject identifier (1-3000).}
#'   \item{time}{Time variable (integer, 0 to 4).}
#'   \item{V}{Time-fixed baseline covariate.}
#'   \item{L}{Time-varying confounder (continuous).}
#'   \item{A}{Time-varying binary exposure.}
#'   \item{M}{Time-varying mediator (continuous).}
#'   \item{Y}{Event indicator at this time point (1 = event, 0 = still at
#'     risk). There is no censoring in this dataset.}
#'   \item{lag1_A}{Exposure at the previous time point; 0 at \code{time == 0}.}
#'   \item{lag1_M}{Previous value of \code{M}; 0 at \code{time == 0}.}
#'   \item{lag1_L}{Previous value of \code{L}; 0 at \code{time == 0}.}
#' }
#'
#' @seealso \code{\link{nonsurvivaldata}} for the end-of-follow-up counterpart.
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
#' The mediated-interaction residual is non-zero because the reported total
#' effect is a natural-course contrast while the direct and indirect effects
#' are interventional; see the \emph{Mediator pool} section of
#' \code{\link{mediation}} for how the mediator draws are constructed.
#'
#' @source Simulated from the data-generating process of Yamamuro, S.,
#'   Shinozaki, T., Iimuro, S., & Matsuyama, Y. (2021). Mediational g-formula
#'   for time-varying treatment and repeated-measured multiple mediators.
#'   \emph{Statistical Methods in Medical Research}, 30(8), 1782-1799.
#'   \doi{10.1177/09622802211025988}
"yamamurodata"
