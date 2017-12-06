##

.LargeLikelihood  <- function(pars, coverage, markerImbalances, ECMCollected, numberOfContributors, logGenotypeProbability, normalisingConstant) {
    numberOfProfiles = dim(ECMCollected)[2] / numberOfContributors

    heterozygoteAverage = pars[1]
    dispersion = pars[2]

    noiseAverage = pars[3]
    noiseDispersion = pars[4]

    phi = rep(pars[5:length(pars)], times = numberOfProfiles)

    markerImbalancesCollected = rep(markerImbalances, times = numberOfProfiles)
    coverageCollected = rep(coverage, times = numberOfProfiles)

    EC = (ECMCollected %*% phi)
    mu_ma = heterozygoteAverage * markerImbalancesCollected * EC

    ll_allele_coverage <- ifelse(EC != 0, dnbinom(coverageCollected, mu = mu_ma, size = dispersion, log = T), 0.0)
    ll_noise <- ifelse(EC == 0, dnbinom(coverageCollected, mu = noiseAverage, size = noiseDispersion, log = T), 0.0)

    ll <- ll_allele_coverage + ll_noise
    llmatrix <- rowSums(matrix(ll, ncol = length(coverage), byrow = T)) + logGenotypeProbability

    LL <- sum(exp(llmatrix + abs(normalisingConstant)))
    return(-LL)
}

# REA = optimalUnknownGenotypesHp; normalisingConstant = likelihoodNormaliser
optimiseParametersLargeLikelihood <- function(sampleTibble, REA, numberOfContributors, normalisingConstant) {
    coverage <- sampleTibble$Coverage
    markerImbalances <- sampleTibble$MarkerImbalance

    ECMFittestProfiles <- lapply(REA, function(rr) rr$ExpectedContributionMatrix)
    ECMLength <- unlist(lapply(ECMFittestProfiles, function(ECM) dim(ECM)[1]))
    ECMLengthPartialSum <- cumsum(c(1, ECMLength))
    ECMCollected <- matrix(0, ncol = numberOfContributors * length(ECMFittestProfiles), nrow = sum(ECMLength))
    for (r in 1:length(ECMFittestProfiles)) {
        ECMCollected[(ECMLengthPartialSum[r]):(ECMLengthPartialSum[r + 1] - 1), ((r - 1) * numberOfContributors + 1):(r * numberOfContributors)] <- ECMFittestProfiles[[r]]
    }

    lowerBounds = c(rep(0, times = 4), rep(0, times = numberOfContributors))
    upperBounds = c(rep(Inf, times = 4), rep(1, times = numberOfContributors))

    eqBounds <- function(pars, coverage, markerImbalances, ECMCollected, numberOfContributors, logGenotypeProbability, normalisingConstant) {
        return(sum(pars[5:length(pars)]))
    }

    logGenotypeProbability = unlist(lapply(REA, function(rr) rr$LogLikelihoods[3]))

    fittestParameters <- REA[[1]]$Parameters
    initialPars <- c(fittestParameters$SampleParameters, fittestParameters$NoiseParameters,
                     fittestParameters$MixtureParameters)

    pars = solnp(initialPars, fun = .LargeLikelihood,
                     LB = lowerBounds, UB = upperBounds, eqfun = eqBounds, eqB = 1,
                     coverage = coverage, markerImbalances = markerImbalances,
                     ECMCollected = ECMCollected, numberOfContributors = numberOfContributors,
                     logGenotypeProbability = logGenotypeProbability, normalisingConstant = normalisingConstant,
                     control = list(trace = FALSE))

    res <- list(SampleParameters = pars$pars[1:2], NoiseParameters = pars$pars[3:4], MixtureParameters = pars$pars[(5):length(pars$pars)], Likelihood = -pars$values[length(pars$values)])
    return(res)
}

