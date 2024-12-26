###############################################################################
###      TEMPORAL TRENDS IN SPATIAL DETERMINANTS OF WOLVERINE DENSITY       ### 
###                          IN NORWAY AND SWEDEN                           ###
###              USING OPEN-POPULATION SPATIAL CAPTURE-RECAPTURE            ###
###                               WITH RJMCMC                               ###
###                                                                         ###
###     Nine monitoring seasons: Dec 2013-Jun 2014 to Dec 2021-Jun 2022     ### 
###                https://doi.org/10.1073/pnas.2401679122                  ###
###############################################################################

## Clean
rm(list = ls())
cat("\014")
gc()

## Library
library(nimble)
library(nimbleSCR) 

## Load a custom nimble function
source("~/dbin_LESS_Cached_MultipleCovResponse.R") 


## -----------------------------------------------------------------------------
## ------ I. LOAD NIMBLE INPUT FILES -----
## -----------------------------------------------------------------------------

# There are two input files, one for female and one for male wolverines.
# The code below assumes data for one sex is loaded only.

## Load wolverine data
load("~/OPSCR_RJ_Female_1.RData") 
load("~/OPSCR_RJ_Male_1.RData")


## -----------------------------------------------------------------------------
## ------ II. MODEL SETTING AND RUNNING ------- 
## -----------------------------------------------------------------------------

### ==== NIMBLE MODEL DEFINITION ====
modelCode <- nimbleCode({

  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##  
  
  # PRIORS
  pi ~ dunif(0,1) #--Prior on inclusion probability
  
  for (h in 1:n.habCovs) {
    beta0[h] ~ dunif(-10,10) #--intercept
    beta1[h] ~ dunif(-10,10) #--slope
    wbeta1[h] ~ dbern(pi)    #--Indicator variable for the temporal effect of each density covariate
  }#h
  
  # REGRESSION MODEL: the intercept and the effect of the year are constant
  for (t in 1:n.years) {
    for (h in 1:n.habCovs) { 
      betaHabCovs[h, t] <- beta0[h] + (wbeta1[h] * beta1[h] * Time[t])
    }#h
    
    # 10 spatial covariate of wolverine density
    habIntensity[1:numHabWindows, t] <- exp(
      
        betaHabCovs[1,  t] * habCovs[1:numHabWindows,  1] + # distRelict 
        betaHabCovs[2,  t] * habCovs[1:numHabWindows,  2] + # ForestCover 
        betaHabCovs[3,  t] * habCovs[1:numHabWindows,  3] + # MooseDensity 
        betaHabCovs[4,  t] * habCovs[1:numHabWindows,  4] + # wildReindeerArea
        betaHabCovs[5,  t] * habCovs[1:numHabWindows,  5] + # sweReindeerArea
        betaHabCovs[6,  t] * habCovs[1:numHabWindows,  6] + # norReindeerArea
        betaHabCovs[7,  t] * habCovs[1:numHabWindows,  7] + # SheepDensity
        betaHabCovs[8,  t] * habCovs[1:numHabWindows,  8] + # Ruggedness
        betaHabCovs[9,  t] * habCovs[1:numHabWindows,  9] + # SnowCovered
        betaHabCovs[10, t] * habCovs[1:numHabWindows, 10]   # Settlements_log
      
    )
    
    sumHabIntensity[t] <- sum(habIntensity[1:numHabWindows, t])
    logHabIntensity[1:numHabWindows, t] <- log(habIntensity[1:numHabWindows, t])
    logSumHabIntensity[t] <- log(sumHabIntensity[t])
    
  }#t
  
  
  for (t in 1:n.years) {
    for (i in 1:n.individuals) {
      sxy[i, 1:2, t] ~ dbernppAC(
        lowerCoords = lowerHabCoords[1:numHabWindows, 1:2]
        ,
        upperCoords = upperHabCoords[1:numHabWindows, 1:2]
        ,
        logIntensities = logHabIntensity[1:numHabWindows, t]
        ,
        logSumIntensity = logSumHabIntensity[t]
        ,
        habitatGrid =  habitatGrid[1:y.max, 1:x.max]
        ,
        numGridRows = y.max
        ,
        numGridCols = x.max
      )
    }#i  
  }#t
  

  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##    
  omeg1[1:3] ~ ddirch(alpha[1:3])   
  
  for (t in 1:n.years1) {
    # PRIORS 
    gamma[t] ~ dunif(0,1)
    phi[t] ~ dunif(0,1)
    
    # "UNBORN"
    omega[1,1,t] <- 1-gamma[t]
    omega[1,2,t] <- gamma[t]
    omega[1,3,t] <- 0
    # "Alive"
    omega[2,1,t] <- 0
    omega[2,2,t] <- phi[t]
    omega[2,3,t] <- 1-phi[t]
    # "Dead"
    omega[3,1,t] <- 0
    omega[3,2,t] <- 0
    omega[3,3,t] <- 1
  }#t
  
  pResponse ~ dunif(0, 1)
  
  for (i in 1:n.individuals) { 
    detResponse[i,1] ~ dbern(pResponse)
    
    z[i,1] ~ dcat(omeg1[1:3]) 
    for (t in 1:n.years1) {
      z[i, t+1] ~ dcat(omega[z[i,t], 1:3, t]) 
    }#i 								
  }#t 
  
  
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------##
  for (t in 1:n.years) {
    sigma[t] ~ dunif(0,4)
    for (c in 1:n.covs) {
      betaCovs[c,t] ~ dunif(-5,5)
    }
    
    for (c in 1:n.covsOth) {
      betaCovsOth[c,t] ~ dunif(-5,5)
    }
    
    betaResponse[t] ~ dunif(-5,5)
    betaResponseOth[t] ~ dunif(-5,5)
  }
  
  for (c in 1:n.counties) {
    for (t in 1:n.years) {
      p01[c,t] ~ dunif(0,1)
      p0[c,t] <- p01[c,t] * countyToggle[c,t]
    }#t
  }#c  
  
  for (c in 1:n.countries) {
    for (t in 1:n.years) {
      p01Oth[c,t] ~ dunif(0,1)
      p0Oth[c,t] <- p01Oth[c,t] * countyToggleOth[c,t]
    }#t
  }#c  
  
  for (t in 1:n.years) {
    for (i in 1:n.individuals) {
      y.alive[i, 1:nMaxDetectors,t] ~ dbin_LESS_Cached_MultipleCovResponse(sxy = sxy[i,1:2,t],
                                                                           sigma = sigma[t],
                                                                           nbDetections[i,t],
                                                                           yDets = yDets[i,1:nMaxDetectors,t],
                                                                           detector.xy =  detector.xy[1:n.detectors,1:2],
                                                                           trials = trials[1:n.detectors],
                                                                           detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
                                                                           nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
                                                                           ResizeFactor = ResizeFactor,
                                                                           maxNBDets = maxNBDets,
                                                                           habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                           indicator = isAlive[i,t],
                                                                           p0[1:n.counties,t],
                                                                           detCounties[1:n.detectors],
                                                                           detCov = detCovs[1:n.detectors,t,1:n.covs],
                                                                           betaCov = betaCovs[1:n.covs,t],
                                                                           BetaResponse = betaResponse[t],
                                                                           detResponse = detResponse[i,t])
      
      y.aliveOth[i, 1:nMaxDetectorsOth,t] ~ dbin_LESS_Cached_MultipleCovResponse(sxy = sxy[i,1:2,t],
                                                                                 sigma = sigma[t],
                                                                                 nbDetectionsOth[i,t],
                                                                                 yDets = yDetsOth[i,1:nMaxDetectorsOth,t],
                                                                                 detector.xy =  detector.xy[1:n.detectors,1:2],
                                                                                 trials = trials[1:n.detectors],
                                                                                 detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
                                                                                 nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
                                                                                 ResizeFactor = ResizeFactor,
                                                                                 maxNBDets = maxNBDets,
                                                                                 habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                                 indicator = isAlive[i,t],
                                                                                 p0Oth[1:n.countries,t],
                                                                                 detCountries[1:n.detectors,t],
                                                                                 detCov = detCovsOth[1:n.detectors,t,1:n.covsOth],
                                                                                 betaCov = betaCovsOth[1:n.covsOth,t],
                                                                                 BetaResponse = betaResponseOth[t],
                                                                                 detResponse = detResponse[i,t])
    }#i
  }#t
  

  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  for (i in 1:n.individuals) { 
    isAlive[i,1] <- (z[i,1] == 2) 
    for (t in 1:n.years1) {
      isAlive[i,t+1] <- (z[i,t+1] == 2) 
    }#t
  }#i
  for (t in 1:n.years) {
    N[t] <- sum(isAlive[1:n.individuals,t])
  }#t
})


### ----------------------------------------------------------------------------
nimParams <- c("N","beta0","beta1","betaHabCovs", 
               "omeg1","gamma","phi","sigma","pResponse","p0","p0Oth",
               "betaCovs","betaResponse","betaCovsOth","betaResponseOth",
               # "habIntensity", 
               "pi","wbeta1") #--RJMCMC related
nimParams2 <- c("z", "sxy")


## -----------------------------------------------------------------------------
## ------ III. MODEL FITTING ------- WITH RJMCMC
## -----------------------------------------------------------------------------

### ==== CONFIGURE NIMBLE MODEL ====
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      inits = nimInits,
                      data = nimData,
                      check = FALSE,
                      calculate = FALSE ) 

# check
system.time(print(model$calculate()))


### ==== RESUME MODEL CONFIGURATION ====
conf <- configureMCMC(model, monitors = nimParams, thin = 1)
configureRJ(conf,
            targetNodes = c('beta1'),
            indicatorNodes = c('wbeta1'),
            control = list(mean = 0, scale = .2))
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = model, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc


### ==== A TEST RUN ====
niters  <- 101
nburns  <- 1
nchains <- 2

nimbleMCMC_samples <- runMCMC( Cmcmc,
                               niter = niters, 
                               nburnin = nburns, 
                               nchains = nchains, 
                               samplesAsCodaMCMC = TRUE )


# traceplot for the estimated trend in the effect of moose density as an example
plot(nimbleMCMC_samples[,'beta1[3]'],  type = 'l', las = 1)
nimbleMCMC_samples[,'wbeta1[3]']
