library(data.table)
dt <- readRDS("data-raw/time_var1.rds")
dt <- data.table(dt)
dt <- dt[, c("lag_A", "lag_LZ") := list(shift(A, fill = 0), shift(LZ, fill = 0)), by = id]
#dt$C <- ifelse(dt$Y == 1, 0, dt$C)
dt <- dt[, cum_Y := Y + shift(Y, fill = 0), by = id]
dt$C <- ifelse(is.na(dt$C), NA,
               ifelse(dt$C == 1, 0, 1))  # Remain uncensored


source("R/iptw.R")
# devtools::load_all()
debug(iptw_med)
res <- iptw_med(data = dt,
                id.vars = "id",
                time.var = "time",
                baseline.vars = c("W1", "W2"),
                exposure.model = "A ~ W2 + lag_A + lag_LZ",
                mediator.model = "Z ~ W2 + A + LA",
                mediator.family = "binomial",
                exposure.vars = c("A", "lag_A"),
                outcome = "Y",
                censor.model = "C ~ W2 + lag_A +lag_LZ",
                bound.value = 0.01,
                R = 10)



