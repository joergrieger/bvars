library(bvars)
data("USMonPol")


#
# Estimate different bayesian VAR models
#

#
# Estimate Bayesian VAR with Minnesota Prior
#

prior  <- set_prior_minnesota(mydata = USMonPol,nolags = 2)
bvestimate <- bvars::bvar(mydata=USMonPol,priorObj = prior,stabletest = T,nthin = 2,nreps = 110, burnin = 10)

#
# Estimate Regime Switching VAR with uninformative prior
#

prior <- set_prior_uninformative(mydata=USMonPol,nolags = 2)
msestimate <- bvars::msvar(mydata=USMonPol,noregimes=2, priorObj = prior,stabletest = T,nthin = 2, nreps = 110, burnin = 10)

#
# Estimate Threshold VAR with ssvs-prior
#

prior <- set_prior_ssvs(mydata=USMonPol,nolags = 1,tau=100,kappa=100)
prior$kappa0 = 0.1
prior$tau0   = 0.1
tvestimate <- tvar(mydata = USMonPol,priorObj = prior,thMax = 10, thVar = 2,nthin = 5, nreps = 110,burnin = 10,stabletest = FALSE)

# Get Impulse-Response-Functions

ident <- set_identification_cholesky()
irfestimate <- bvars::irf(bvestimate, id_obj = ident, nhor = 24, ncores = 1)

msirfestimate <- bvars::irf(msestimate,id_obj=ident,nhor=12,ncores=1)
tvirestimate  <- bvars::irf(tvestimate,id_obj=ident,nhor=12,ncores=2,bootrep=1)


#
# Plot Impulse-Response Functions
#

plot(irfestimate)
plot(msirfestimate)
plot(tvirestimate)

#
# test forecasts
#

fc <- forecast(bvestimate,forecastHorizon = 6)
plot(fc)

fctv <- forecast(tvestimate,forecastHorizon = 6)
plot(fctv)

msfc <- forecast(msestimate,forecastHorizon = 6)
plot(msfc)

#
# test historical decomposition
#
#hd_1 <- hd(bvestimate)
#hd_2 <- hd(msestimate)
hd_3 <- hd(tvestimate)



#
# Test forecast error variance decomposition
#

fevd_1 <- fevd(bvestimate)
fevd_2 <- fevd(msestimate)
fevd_3 <- fevd(tvestimate)
#

