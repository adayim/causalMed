dat <- readRDS("data/data.rds")
dat$status <- as.numeric(dat$status)
dt <- dat[dat$tij == 0, ]

res <- iorw(data = dt, trt = "a", med = "mt",
            y = "status", time = "eventtime",
            family = "cox",
            cov = c("w1", "w2"))

res <- medlong(data = dat,
               trt = "a",
               med = "mt",
               y = "status",
               id = "id",
               time = "tij",
               cov = c("w1", "w2"),
               m.family  = "gaussian",
               y.family  = "binomial")


