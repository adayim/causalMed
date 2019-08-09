load("data-raw/data_raw.RData")

set.seed(20190809)

dummy <- dummy[dummy$age == dummy$age0, ]
dummy <- dummy[grep("1003_", dummy$pcode), ]
dummy$id <- gsub("1003_|ID", "1", dummy$pcode)
dummy$smoke <- dummy$imp_smoke

tmp <- aggregate(age ~ pcode, data = dat, min)
tmp$age0 <- tmp$age
tmp$age <- NULL

dat <- merge(tmp, dat, by = "pcode")
dat$id <- gsub("1003_|ID", "1", dat$pcode)
dat <- dat[, c("id", "age0", "hdl", "ldl", "tg", "age")]

dfm <- dummy[, c("id", "age0", "gender", "class", "smoke", "bmi", "os", "cvd")]

dfm <- merge(dfm, dat, by = c("id", "age0"), all = TRUE)
dfm$time <- dfm$age - dfm$age0
lipdat <- transform(dfm, id=as.numeric(factor(id)))
lipdat$age <- NULL
lipdat <- lipdat[, c("id", "age0", "gender",
                     "class", "smoke", "bmi",
                     "hdl", "ldl", "tg", "time",
                     "cvd", "os")]
lipdat <- lipdat[lipdat$id %in% sample(lipdat$id, 500), ]

devtools::use_data(lipdat, overwrite = TRUE)


