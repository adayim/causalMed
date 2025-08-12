#' Lipid test results of 500 person.
#'
#' A dataset containing the lipid profile and other attributes.
#'
#' @format A data frame with 2114 rows and 13 variables:
#' \describe{
#'   \item{id}{unique id of subject}
#'   \item{age}{baseline age}
#'   \item{gender}{gender of subject}
#'   \item{class}{class membership of per-subject}
#'   \item{smoke}{Ever smoke, 0=No, 1=Yes}
#'   \item{bmi}{BMI, in kg/^2}
#'   \item{hdl}{high density lipoprotein, mmol/L, time-varying}
#'   \item{ldl}{low density lipoprotein, mmol/L, time-varying}
#'   \item{tg}{triglyceride, mmol/L, time-varying}
#'   \item{time}{Followup time since baseline(time = 0)}
#'   \item{cvd}{cardiovascular diseases event}
#'   \item{os}{Overall Survival of cardiovascular diseases}
#' }
#' @source \url{http://www.diamondse.info/}
"lipdat"


#' Bone marrow transplant data
#'
#' The study population arose from a multicenter trial of leukemia patients and
#' comprises 137 individuals prepared for bone marrow transplants under a radiation-free regimen at four medical centers.
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
#'   \item{d}{indicator of death (1= yes, 0=no) }
#'   \item{gvhd}{indicator of GvHD (1= yes, 0=no)}
#'   \item{daysgvhd}{number of days since onset of GvHD}
#'   \item{daysnogvhd}{number of days without GvHD}
#'   \item{gvhdm1}{lagged GvHD}
#'   \item{relapse}{relapse indicator}
#'   \item{relapsem1}{lagged relapse indicator}
#'   \item{daysrelapse}{number of days since relapse}
#'   \item{daysnorelapse}{number of days relapse-free}
#'   \item{platnorm}{normal platelet levels (1=patient  has relapsed or reached
#'   normal platelet count, 0=not in relapse or below normal  platelets)}
#'   \item{platnormm1}{lagged normal platelet levels}
#'   \item{daysplatnorm}{time spent reaching normal platelet levels}
#'   \item{daysnoplatnorm}{time spent without reaching normal platelet levels}
#'   \item{censlost}{indicator of censoring due to loss-to-follow up (1=yes, 0=no)}
#' }
#' @references Keil et al. (2014) Epidemiology Epidemiology, 25(6), 889–897
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4310506/}{PubMed})
"gvhd"
