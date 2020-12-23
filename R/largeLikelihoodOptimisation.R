##
.itpg <- function(x, mu, theta, p, t, i) {
    .tpg <- function(x, mu, theta, t) {
        return(dnbinom(x, mu = mu, size = theta, log = T) - log(1 - pnbinom(t, mu = mu, size = theta, log = F)))
    }

    x_i <- (x == i)
    xx <- rep(0, length(x))
    xx[x_i] <- log(p + exp(log(1.0 - p) + .tpg(i, mu, theta, t)))
    xx[!x_i] <- log(1.0 - p) + .tpg(x[!x_i], mu, theta, t)

    return(xx)
}

.LargeLikelihood  <- function(pars, coverage, markerImbalances, ECMCollected, numberOfContributors, logGenotypeProbability, normalisingConstant) {
    numberOfProfiles = length(ECMCollected)

    heterozygoteAverage = pars[1]
    dispersion = pars[2]

    noiseAverage = pars[3]
    noiseDispersion = pars[4]
    noiseInflation = pars[5]

    phi = pars[6:length(pars)]

    ll <- rep(0, numberOfProfiles)
    for (i in 1:numberOfProfiles) {
        EC = (ECMCollected[[i]] %*% phi)
        mu_ma = heterozygoteAverage * markerImbalances * EC

        alleles <- EC != 0

        ll_allele_coverage <- dnbinom(coverage[alleles], mu = mu_ma[alleles], size = mu_ma[alleles] / dispersion, log = T)
        ll_noise <- .itpg(coverage[!alleles], mu = noiseAverage, theta = noiseDispersion, p = noiseInflation, t = 0.0, i = 1.0)

        ll[i] <- sum(ll_allele_coverage) + sum(ll_noise) + logGenotypeProbability[i]
    }


    LL <- sum(exp(ll + abs(normalisingConstant)))
    return(-LL)
}

# REA = optimalUnknownGenotypesHp; normalisingConstant = likelihoodNormaliser
.optimiseParametersLargeLikelihood <- function(sampleTibble, REA, numberOfContributors, normalisingConstant) {
    coverage <- sampleTibble$Coverage

    ECMCollected <- lapply(REA, function(rr) rr$ExpectedContributionMatrix)

    lowerBounds = c(rep(0, times = 5), rep(0, times = numberOfContributors))
    lowerBounds[4] = 1.0

    upperBounds = c(rep(Inf, times = 5), rep(1, times = numberOfContributors))
    upperBounds[5] = 1.0

    eqBounds <- function(pars, coverage, markerImbalances, ECMCollected, numberOfContributors, logGenotypeProbability, normalisingConstant) {
        return(sum(pars[6:length(pars)]))
    }

    logGenotypeProbability = unlist(lapply(REA, function(rr) rr$LogLikelihoods[3]))

    fittestParameters <- REA[[1]]$Parameters
    initialPars <- c(fittestParameters$SampleParameters, fittestParameters$NoiseParameters,
                     fittestParameters$MixtureParameters)

    markerImbalances = (
        sampleTibble %>% left_join(sampleTibble %>%
                                       distinct(Marker) %>%
                                       mutate(MI = fittestParameters$MarkerImbalanceParameters), by = "Marker")
    )$MI

    pars = solnp(initialPars, fun = .LargeLikelihood,
                 LB = lowerBounds, UB = upperBounds, eqfun = eqBounds, eqB = 1,
                 coverage = coverage, markerImbalances = markerImbalances,
                 ECMCollected = ECMCollected, numberOfContributors = numberOfContributors,
                 logGenotypeProbability = logGenotypeProbability, normalisingConstant = normalisingConstant,
                 control = list(trace = FALSE))

    res <- list(SampleParameters = pars$pars[1:2], NoiseParameters = pars$pars[3:5], MixtureParameters = pars$pars[6:length(pars$pars)],
                MarkerImbalanceParameters = fittestParameters$MarkerImbalanceParameters, Likelihood = -pars$values[length(pars$values)])
    return(res)
}

