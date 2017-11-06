#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//[[Rcpp::depends(BH)]]
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include "logLikelihoods.hpp"
#include "AuxiliaryFunctions.hpp"

double logPoissonGammaDistribution(double & x, double & mean, double & dispersion)
{
    double logProbability = boost::math::lgamma(x + dispersion) - boost::math::lgamma(x + 1) - boost::math::lgamma(dispersion) +
        x * (std::log(mean) - std::log(mean + dispersion)) + dispersion * (std::log(dispersion) - std::log(mean + dispersion));
    return logProbability;
}

double logMultinomialCoefficient(int totalCounts, Eigen::VectorXd counts)
{
    std::size_t N = counts.size();

    double logGammaCounts = 0;
    for (std::size_t n = 0; n < N; n++)
    {
        logGammaCounts += lgamma(counts[n] + 1);
    }
    return lgamma(totalCounts + 1) - logGammaCounts;
}

double logDirichletMultinomial(Eigen::VectorXd counts, Eigen::VectorXd alleleFrequencies)
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

double logDirichletMultinomialTheta(Eigen::VectorXd counts, double theta, Eigen::VectorXd alleleFrequencies)
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

double logLikelihoodAlleleCoverage(Eigen::VectorXd coverage, Eigen::MatrixXd expectedContributionMatrix,
                                   Eigen::VectorXd sampleParameters, Eigen::VectorXd mixtureProportions,
                                   Eigen::VectorXd markerImbalances)
{
    Eigen::VectorXd expectedContribution = expectedContributionMatrix * mixtureProportions;

    double referenceMarkerAverage = sampleParameters[0];
    double dispersion = sampleParameters[1];

    double logLikelihood = 0.0;
    for (std::size_t n = 0; n < coverage.size(); n++)
    {
        double mu_ma = referenceMarkerAverage * markerImbalances[n] * expectedContribution[n];
        if (mu_ma > 0.0)
        {
            logLikelihood += logPoissonGammaDistribution(coverage[n], mu_ma, dispersion);
        }
    }

    return logLikelihood;
}

double logLikelihoodNoiseCoverage(Eigen::VectorXd & coverage, Eigen::VectorXd & noiseProfile, double noiseLevel, double noiseDispersion)
{
    std::size_t N = coverage.size();

    double logLikelihood = 0.0;
    for (std::size_t n = 0; n < N; n++)
    {
        if (noiseProfile[n] == 1.0)
        {
            logLikelihood += logPoissonGammaDistribution(coverage[n], noiseLevel, noiseDispersion);
        }
    }

    return logLikelihood;
}

double logGenotypeProbabilityHWE(Eigen::VectorXd & alleleFrequencies, Eigen::MatrixXd & unknownProfiles, std::size_t numberOfMarkers, Eigen::VectorXd & numberOfAlleles)
{
    Eigen::VectorXd partialNumberOfAlleles = partialSumEigen(numberOfAlleles);

    double logPriorProbability = 0.0;
    for (std::size_t m = 0; m < numberOfMarkers; m++)
    {
        Eigen::MatrixXd profiles_m = unknownProfiles.block(partialNumberOfAlleles[m], 0, numberOfAlleles[m], unknownProfiles.cols());

        Eigen::VectorXd profileCounts_m = profiles_m.rowwise().sum();
        Eigen::VectorXd alleleFractions_m = alleleFrequencies.segment(partialNumberOfAlleles[m], numberOfAlleles[m]);

        logPriorProbability += logDirichletMultinomial(profileCounts_m, alleleFractions_m);
    }

    return logPriorProbability;
}

double logGenotypeProbabilityThetaCorrection(Eigen::VectorXd & alleleFrequencies, double theta, Eigen::MatrixXd & unknownProfiles, Eigen::MatrixXd & knownProfiles, std::size_t numberOfMarkers, Eigen::VectorXd numberOfAlleles)
{
    Eigen::MatrixXd profiles = bindColumns(knownProfiles, unknownProfiles);
    std::size_t numberOfContributors = profiles.cols();
    std::size_t numberOfKnownContributors = numberOfContributors - unknownProfiles.cols();

    Eigen::VectorXd partialNumberOfAlleles = partialSumEigen(numberOfAlleles);
    double thetaCorrectedProbabilty = 0.0;
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

        thetaCorrectedProbabilty += logProbabilty_m;
    }

    return thetaCorrectedProbabilty;
}

//[[Rcpp::export()]]
double logPriorGenotypeProbability(Eigen::VectorXd & alleleFrequencies, double theta, Eigen::MatrixXd & unknownProfiles, Eigen::MatrixXd & knownProfiles, std::size_t numberOfMarkers, Eigen::VectorXd numberOfAlleles)
{
    double logPriorProbability = 0.0;
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
