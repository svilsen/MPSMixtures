#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//[[Rcpp::depends(BH)]]
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include "logLikelihoods.hpp"
#include "AuxiliaryFunctions.hpp"

double logPoissonGammaDistribution(const double & x, const double & mean, const double & dispersion)
{
    double logProbability = boost::math::lgamma(x + dispersion) - boost::math::lgamma(x + 1) - boost::math::lgamma(dispersion) +
        x * (std::log(mean) - std::log(mean + dispersion)) + dispersion * (std::log(dispersion) - std::log(mean + dispersion));
    return logProbability;
}

//[[Rcpp::export(.devianceResidualPoissonGammaDistribution)]]
double devianceResidualPoissonGammaDistribution(const double & x, const double & mean, const double & dispersion)
{
	if (std::abs(x - mean) < 2e-8)
	{
		return 0.0;
	}

    double deviance_ma = (x + dispersion) * (std::log(mean + dispersion) - std::log(x + dispersion));
    if (x > 0)
    {
        deviance_ma += x * (std::log(x) - std::log(mean));
    }

    return boost::math::sign(x - mean) * std::pow(2.0 * deviance_ma, 0.5);
}

double devianceResidualPG1(const double & x, const double & mean, const double & dispersion)
{
    double deviance_ma = (x - mean) * std::log(dispersion + 1) / dispersion - boost::math::lgamma((dispersion * x + mean) / dispersion) +
        boost::math::lgamma(x * (dispersion + 1) / dispersion);
    if (x > 0)
    {
        deviance_ma += boost::math::lgamma(x / dispersion) - boost::math::lgamma(mean / dispersion);
    }

    return boost::math::sign(x - mean) * std::pow(2.0 * deviance_ma, 0.5);
}


double logMultinomialCoefficient(const int & totalCounts, const Eigen::VectorXd & counts)
{
    std::size_t N = counts.size();

    double logGammaCounts = 0;
    for (std::size_t n = 0; n < N; n++)
    {
        logGammaCounts += lgamma(counts[n] + 1);
    }
    return lgamma(totalCounts + 1) - logGammaCounts;
}

double logDirichletMultinomial(const Eigen::VectorXd & counts, const Eigen::VectorXd & alleleFrequencies)
{
    std::size_t N = alleleFrequencies.size();
    int sumCounts = counts.sum();

    double groupSum = 0.0;
    for (std::size_t n = 0; n < N; n++)
    {
        groupSum += counts[n] * std::log(alleleFrequencies[n]);

    }

    double resLogDirichletMultinomial = logMultinomialCoefficient(sumCounts, counts) + groupSum;
    return resLogDirichletMultinomial;
}

double logDirichletMultinomialTheta(const Eigen::VectorXd & counts, const double & theta, const Eigen::VectorXd & alleleFrequencies)
{
    std::size_t N = alleleFrequencies.size();
    int sumCounts = counts.sum();

    double groupSum = 0.0;
    double gamma_plus = 1/theta - 1;
    for (std::size_t n = 0; n < N; n++)
    {
        double gamma_j = gamma_plus * alleleFrequencies[n];
        groupSum += lgamma(counts[n] + gamma_j) - lgamma(gamma_j);
    }

    groupSum += (lgamma(gamma_plus) - lgamma(sumCounts + gamma_plus));

    double resLogDirichletMultinomial = logMultinomialCoefficient(sumCounts, counts) + groupSum;
    return resLogDirichletMultinomial;
}

Eigen::VectorXd logLikelihoodAlleleCoverage(const Eigen::VectorXd & coverage, const std::vector<Eigen::MatrixXd> & expectedContributionMatrix,
                                            const std::vector<Eigen::VectorXd> & alleleIndex, const Eigen::VectorXd partialSumAlleles,
                                            const Eigen::VectorXd & sampleParameters, const Eigen::VectorXd & mixtureProportions,
                                            const Eigen::VectorXd & markerImbalances)
{
    const std::size_t & M = alleleIndex.size();

    const double & referenceMarkerAverage = sampleParameters[0];
    const double & dispersion = sampleParameters[1];

    Eigen::VectorXd logLikelihood = Eigen::VectorXd::Zero(M);
    for (std::size_t m = 0; m < M; m++)
    {
        const Eigen::VectorXd & alleleIndex_m = alleleIndex[m];
        const Eigen::MatrixXd & expectedContributionMatrix_m = expectedContributionMatrix[m];
        const Eigen::VectorXd & expectedContribution = expectedContributionMatrix_m * mixtureProportions;

        double logLikelihood_m = 0.0;
        for (std::size_t a = 0; a < alleleIndex_m.size(); a++)
        {
            const std::size_t & n = partialSumAlleles[m] + alleleIndex_m[a];

            double mu_ma = referenceMarkerAverage * markerImbalances[m] * expectedContribution[a];
            if (expectedContribution[a] > 0)
            {
                logLikelihood_m += logPoissonGammaDistribution(coverage[n], mu_ma, mu_ma / dispersion);
            }
        }

        logLikelihood[m] = logLikelihood_m;
    }

    return logLikelihood;
}

Eigen::VectorXd logLikelihoodNoiseCoverage(const Eigen::VectorXd & coverage, const std::vector<Eigen::VectorXd> & noiseIndex,
                                           const Eigen::VectorXd & partialSumAlleles, const double & noiseLevel, const double & noiseDispersion)
{
    const std::size_t & M = noiseIndex.size();
    const double & dispersion = noiseDispersion;

    Eigen::VectorXd logLikelihood = Eigen::VectorXd::Zero(M);
    for (std::size_t m = 0; m < M; m++)
    {
        const Eigen::VectorXd & noiseIndex_m = noiseIndex[m];
        for (std::size_t a = 0; a < noiseIndex_m.size(); a++)
        {
            const std::size_t & n = partialSumAlleles[m] + noiseIndex_m[a];
            logLikelihood[m] += logPoissonGammaDistribution(coverage[n], noiseLevel, dispersion) -
                std::log(1 - std::exp(dispersion * (std::log(dispersion) - std::log(noiseLevel + dispersion))));
        }
    }

    return logLikelihood;
}


Eigen::VectorXd logGenotypeProbabilityHWE(const Eigen::VectorXd & alleleFrequencies, const Eigen::MatrixXd & unknownProfiles, const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles)
{
    Eigen::VectorXd partialNumberOfAlleles = partialSumEigen(numberOfAlleles);

    Eigen::VectorXd logPriorProbability = Eigen::VectorXd::Zero(numberOfMarkers);
    for (std::size_t m = 0; m < numberOfMarkers; m++)
    {
        Eigen::MatrixXd profiles_m = unknownProfiles.block(partialNumberOfAlleles[m], 0, numberOfAlleles[m], unknownProfiles.cols());

        Eigen::VectorXd profileCounts_m = profiles_m.rowwise().sum();
        Eigen::VectorXd alleleFractions_m = alleleFrequencies.segment(partialNumberOfAlleles[m], numberOfAlleles[m]);

        logPriorProbability[m] = logDirichletMultinomial(profileCounts_m, alleleFractions_m);
    }

    return logPriorProbability;
}


Eigen::VectorXd logGenotypeProbabilityThetaCorrection(const Eigen::VectorXd & alleleFrequencies, const double & theta, const Eigen::MatrixXd & unknownProfiles, const Eigen::MatrixXd & knownProfiles, const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles)
{
    Eigen::MatrixXd profiles = bindColumns(knownProfiles, unknownProfiles);
    std::size_t numberOfContributors = profiles.cols();
    std::size_t numberOfKnownContributors = numberOfContributors - unknownProfiles.cols();

    Eigen::VectorXd partialNumberOfAlleles = partialSumEigen(numberOfAlleles);
    Eigen::VectorXd thetaCorrectedProbabilty = Eigen::VectorXd::Zero(numberOfMarkers);
    for (std::size_t m = 0; m < numberOfMarkers; m++)
    {
        Eigen::MatrixXd profiles_m = profiles.block(partialNumberOfAlleles[m], 0, numberOfAlleles[m], numberOfContributors);

        Eigen::VectorXd profileCounts = profiles_m.rowwise().sum();
        Eigen::VectorXd alleleFractions_m = alleleFrequencies.segment(partialNumberOfAlleles[m], numberOfAlleles[m]);

        double logProbabilityOfUnknownProfiles = logDirichletMultinomialTheta(profileCounts, theta, alleleFractions_m);
        double logProbabilty_m = logProbabilityOfUnknownProfiles;

        if (numberOfKnownContributors != 0)
        {
            Eigen::MatrixXd knownProfiles_m = knownProfiles.block(partialNumberOfAlleles[m], 0, numberOfAlleles[m], numberOfKnownContributors);

            Eigen::VectorXd knownCounts = knownProfiles_m.rowwise().sum();

            double logProbabilityOfKnownProfiles = logDirichletMultinomialTheta(knownCounts, theta, alleleFractions_m);
            logProbabilty_m -= logProbabilityOfKnownProfiles;
        }

        thetaCorrectedProbabilty[m] = logProbabilty_m;
    }

    return thetaCorrectedProbabilty;
}

Eigen::VectorXd logPriorGenotypeProbability(const Eigen::VectorXd & alleleFrequencies, const double & theta, const Eigen::MatrixXd & unknownProfiles, const Eigen::MatrixXd & knownProfiles, const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles)
{
    Eigen::VectorXd logPriorProbability;
    if (theta == 0.0)
    {
        logPriorProbability = logGenotypeProbabilityHWE(alleleFrequencies, unknownProfiles, numberOfMarkers, numberOfAlleles);
    }
    else
    {
        logPriorProbability = logGenotypeProbabilityThetaCorrection(alleleFrequencies, theta, unknownProfiles, knownProfiles, numberOfMarkers, numberOfAlleles);
    }

    return logPriorProbability;
}
