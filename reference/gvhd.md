# Bone marrow transplant data (graft-versus-host disease)

Person-day data for the 137 leukemia patients who received bone marrow
transplants under a radiation-free regimen at four medical centers, as
used in the parametric g-formula illustration of Keil et al. (2014).
Each subject contributes one row per day of follow-up (day 1 up to
death, loss to follow-up, or day 1825), giving 108,714 rows. The dataset
reproduces the person-day file created in Appendix 1 of that paper: 80
deaths and 39 losses to follow-up.

## Usage

``` r
gvhd
```

## Format

A data frame with 108,714 rows and 21 variables:

- id:

  unique subject identifier (1-137)

- age:

  baseline age in years

- male:

  male indicator (1 = male, 0 = female)

- cmv:

  baseline cytomegalovirus immune status (1 = positive, 0 = negative)

- all:

  leukemia type (1 = acute lymphocytic, 0 = acute myeloid)

- wait:

  waiting time from leukemia diagnosis to transplantation, in months
  (waiting days / 30.5)

- day:

  days since bone marrow transplant (1-1825)

- d:

  indicator of death at the end of day `day` (1 = yes, 0 = no)

- gvhd:

  indicator of graft-versus-host disease having occurred by day `day` (1
  = yes, 0 = no); the exposure

- daysgvhd:

  cumulative number of days with GvHD

- daysnogvhd:

  cumulative number of days without GvHD

- gvhdm1:

  one-day lag of `gvhd`

- relapse:

  indicator of leukemia relapse having occurred by day `day` (1 = yes, 0
  = no)

- relapsem1:

  one-day lag of `relapse`

- daysrelapse:

  cumulative number of days since relapse

- daysnorelapse:

  cumulative number of days relapse-free

- platnorm:

  indicator of platelet recovery having occurred by day `day` (1 =
  platelets returned to normal levels, 0 = not yet)

- platnormm1:

  one-day lag of `platnorm`

- daysplatnorm:

  cumulative number of days with normal platelet levels

- daysnoplatnorm:

  cumulative number of days without normal platelet levels

- censlost:

  indicator of censoring due to loss to follow-up at day `day` (1 = yes,
  0 = no)

## Details

GvHD (`gvhd`) is the exposure, death (`d`) is the outcome, and relapse
(`relapse`) and platelet recovery (`platnorm`) are time-varying
confounders. The three time-varying states are *absorbing* (once 1, they
remain 1); their `daysX`/`daysnoX` counters and one-day lags (`Xm1`)
summarize the history entering each model.

A worked total-effect analysis on this dataset (counterfactual mortality
risk had GvHD never occurred versus the natural course), mirroring the
model specifications in Appendix 2 of Keil et al. (2014), is walked
through in
[`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md).

## References

Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., & Cole, S.
R. (2014). The parametric g-formula for time-to-event data: intuition
and a worked example. *Epidemiology*, 25(6), 889-897.
[doi:10.1097/EDE.0000000000000160](https://doi.org/10.1097/EDE.0000000000000160)
