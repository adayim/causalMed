
library(devtools)
load_all()
data(nonsurvivaldata)
print("Unique time values:")
ts <- unique(nonsurvivaldata$time)
print(ts)
print(paste("Length:", length(ts)))
