
library(devtools)
load_all()

data(nonsurvivaldata)
print("Names of nonsurvivaldata:")
print(names(nonsurvivaldata))

print("Is 'Y' in nonsurvivaldata?")
print("Y" %in% names(nonsurvivaldata))

# Simulate what gformula does
models <- list(
  spec_model(Y ~ A + L + V, var_type="binomial", mod_type="outcome")
)
data <- nonsurvivaldata
id_var = "id"
base_vars = "V"
exposure = "A"
time_var = "time"
outcome = "Y"

print("Checking check_var_in:")
tryCatch({
  causalMed:::check_var_in(c(id_var, base_vars, exposure, outcome, time_var), data)
  print("check_var_in passed")
}, error = function(e) {
  print(paste("check_var_in failed:", e$message))
})
