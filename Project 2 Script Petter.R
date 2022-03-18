library("INLA")
load(file='rain.rda') # loads 'rain' data.frame

?control.fixed
?f
inla.doc("rw1")

## a)

hyper = list ( prec = list (prior = "gaussian",
                            param = ?)
control.inla = list(strategy="simplified.laplace", int.strategy="ccd")

ptm <- proc.time()

mod <- inla(n.rain ~ -1 + f(day, model="rw1", constr=FALSE, hyper = hyper),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)

proc.time()[3] - ptm[3]

mod$summary.fitted.values

## b)

## c)

# has intercept

mod2 <- inla(n.rain ~ f(day, model="rw1", constr=TRUE),
            data=rain, Ntrials=n.years, control.compute=list(config = TRUE),
            family="binomial", verbose=TRUE, control.inla=control.inla)
