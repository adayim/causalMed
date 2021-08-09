# normal
dat <- foreign::read.dta("https://stats.idre.ucla.edu/stat/data/truncreg.dta")

mod_cov1 <- spec_model(achiv ~ langscore + prog,
                       var_type = "normal",
                       mod_type = "covriate")

mod_cov1$call$data <- substitute(dat)
fit <- eval(mod_cov1$call)
summary(fit)



# Custom function, eg: truncreg
dat <- foreign::read.dta("https://stats.idre.ucla.edu/stat/data/truncreg.dta")
m <- truncreg::truncreg(achiv ~ langscore + prog, data = dat,
                        point = 40, direction = "left")

mod_cov1 <- spec_model(achiv ~ langscore + prog,
                       family = "truncreg",
                       type = "covriate",
                       custom_fn = truncreg::truncreg,
                       point = 40, direction = "left")

mod_cov1$call$data <- substitute(dat)
fit <- eval(mod_cov1$call)
summary(fit)

# Multinomial
library(nnet)
ml <- foreign::read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")
ml$prog2 <- relevel(ml$prog, ref = "academic")
fit <- multinom(prog2 ~ ses + write, data = ml)
pred <- predict(fit, newdata = ml, type =  "probs")
hm <- Hmisc::rMultinom(pred, 1)
rm <- rmultinom(1, 1, pred)

mod_cov1 <- spec_model(prog2 ~ ses + write,
                       family = "multinomial",
                       type = "covriate",
                       #subset = platnormm1 == 0
)

mod_cov1$call$data <- substitute(ml)
fit <- eval(mod_cov1$call)
summary(fit)


