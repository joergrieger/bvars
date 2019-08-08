library(bvar2)
data("USMonPol")


#
# Estimate different bayesian VAR models
#

#
# Estimate Bayesian VAR with Minnesota Prior
#

prior  <- set_prior_minnesota(mydata = USMonPol,nolags = 2)
bvestimate <- bvar2::bvar(mydata=USMonPol,priorObj = prior,stabletest = T,nthin = 2,nreps = 1100, burnin = 100)

#
# Estimate Regime Switching VAR with uninformative prior
#

prior <- set_prior_uninformative(mydata=USMonPol,nolags = 2)
msestimate <- msvar(mydata=USMonPol,noregimes=2, priorObj = prior,stabletest = T,nthin = 2, nreps = 1100, burnin = 100)

#
# Estimate Threshold VAR with ssvs-prior
#

prior <- set_prior_ssvs(mydata=USMonPol,nolags=3,tau=100,kappa=100)
prior$kappa0 = 0.1
prior$tau0   = 0.1
tvestimate <- tvar(mydata = USMonPol,priorObj = prior,thMax = 10, thVar = 2,nthin = 5, nreps = 110,burnin = 10)

# Get Impulse-Response-Functions

ident <- set_identification_cholesky()
irfestimate <- bvar2::irf(bvestimate, id_obj = ident, nhor = 24, ncores = 2)

msirfestimate <- bvar2::irf(msestimate,id_obj=ident,nhor=12,ncores=1)

tvirestimate  <- bvar2::irf(tvestimate,id_obj=ident, nhor=12,ncores = 1, bootrep = 1)
#
# Plot Impulse-Response Functions
#

plot(irfestimate)
plot(msirfestimate)


#
# test forecasts
#

fc <- forecast(bvestimate,forecastHorizon = 6)
plot(fc)

fctv <- forecast(tvestimate,forecastHorizon = 6)
plot(fctv)

msfc <- forecast(msestimate,forecastHorizon = 6)
plot(msfc)
