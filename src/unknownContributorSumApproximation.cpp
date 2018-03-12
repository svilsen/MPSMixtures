#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include "experimentalSetup.hpp"
#include "individual.hpp"
#include "encodingScheme.hpp"
#include "logLikelihoods.hpp"

#include "AuxiliaryFunctions.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

//[[Rcpp::export(.oneStepApproximationCpp)]]
Rcpp::List oneStepApproximationCpp(const Eigen::VectorXd & encodedProfiles, const Eigen::VectorXd & sampleParameters,
                                   const Eigen::VectorXd & noiseParameters, const Eigen::VectorXd & mixtureParameters,
                                   const Eigen::VectorXd & coverage, const Eigen::VectorXd markerImbalances,
                                   const std::vector< std::vector< Eigen::MatrixXd > > potentialParents,
                                   const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                                   const Eigen::VectorXd & alleleFrequencies, const double & theta,
                                   const std::size_t & numberOfContributors, const std::size_t & numberOfMarkers,
                                   const Eigen::VectorXd & numberOfAlleles,
                                   const std::size_t & levelsOfStutterRecursion)
{
    const std::size_t numberOfKnownContributors = knownProfiles.cols();
    const ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors,
                               knownProfiles, allKnownProfiles, coverage, potentialParents, markerImbalances, 0.8,
                               1e-6, theta, alleleFrequencies, levelsOfStutterRecursion);

    const std::size_t & N = numberOfAlleles[0] * (numberOfAlleles[0] + 1) / 2.0;
    const std::size_t & numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;

    std::vector<Eigen::MatrixXd> genotypeMatrices(N);
    Eigen::VectorXd estimatedProbability = Eigen::VectorXd::Zero(N);
    double estimatedProbabilitySum = 0.0;
    double maxValue = -HUGE_VAL;

    std::size_t n = 0;
    for (std::size_t u = 0; u < numberOfUnknownContributors; u++)
    {
        const std::size_t & a = 2 * u;
        Eigen::VectorXd encodedProfiles_um = encodedProfiles;
        for (std::size_t i = 0; i < numberOfAlleles[0]; i++)
        {
            Eigen::VectorXd encodedProfiles_um = encodedProfiles;
            double tt = encodedProfiles_um[a] + 1;
            encodedProfiles_um[a] = i;

            for (std::size_t j = i; j < numberOfAlleles[0]; j++)
            {

                encodedProfiles_um[a + 1] = j;
                Individual I_umij(encodedProfiles_um, sampleParameters, noiseParameters, mixtureParameters, markerImbalances, ES);

                genotypeMatrices[n] = decoding(I_umij.EncodedProfile, ES.NumberOfAlleles, 1, numberOfUnknownContributors);
                estimatedProbability[n] = I_umij.Fitness;
                estimatedProbabilitySum += std::exp(I_umij.Fitness);

                if (I_umij.Fitness > maxValue)
                    maxValue = I_umij.Fitness;

                n++;
            }
        }
    }

    Eigen::VectorXd normalisedProbability = (estimatedProbability - maxValue * Eigen::VectorXd::Ones(N)).array().exp();
    normalisedProbability = normalisedProbability / normalisedProbability.sum();

    return Rcpp::List::create(Rcpp::Named("GenotypeMatrix") = genotypeMatrices,
                              Rcpp::Named("LogUnnormalisedProbabilitySum") = std::log(estimatedProbabilitySum),
                              Rcpp::Named("LogUnnormalisedProbability") = estimatedProbability,
                              Rcpp::Named("NormalisedProbabilities") = normalisedProbability);
}


//[[Rcpp::export(.EAApproximationCpp)]]
Rcpp::List EAApproximationCpp(const std::vector<Eigen::VectorXd> & encodedProfiles, const Eigen::VectorXd & sampleParameters,
                              const Eigen::VectorXd & noiseParameters, const Eigen::VectorXd & mixtureParameters,
                              const Eigen::VectorXd & coverage, const Eigen::VectorXd markerImbalances,
                              const std::vector< std::vector< Eigen::MatrixXd > > potentialParents,
                              const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                              const Eigen::VectorXd & alleleFrequencies, const double & theta,
                              const std::size_t & numberOfContributors, const std::size_t & numberOfMarkers,
                              const Eigen::VectorXd & numberOfAlleles,
                              const std::size_t & levelsOfStutterRecursion)
{
    const std::size_t numberOfKnownContributors = knownProfiles.cols();
    const ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors,
                               knownProfiles, allKnownProfiles, coverage, potentialParents, markerImbalances, 0.8,
                               1e-6, theta, alleleFrequencies, levelsOfStutterRecursion);

    const std::size_t & numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;
    const std::size_t & N = encodedProfiles.size();

    std::vector<Eigen::MatrixXd> genotypeMatrices(N);
    Eigen::VectorXd estimatedProbability = Eigen::VectorXd::Zero(N);
    double estimatedProbabilitySum = 0.0;
    double maxValue = -HUGE_VAL;
    for (std::size_t n = 0; n < N; n++)
    {
        const Eigen::VectorXd & encodedProfiles_n = encodedProfiles[n];
        Individual I(encodedProfiles_n, sampleParameters, noiseParameters, mixtureParameters, markerImbalances, ES);

        genotypeMatrices[n] = decoding(I.EncodedProfile, ES.NumberOfAlleles, 1, numberOfUnknownContributors);
        estimatedProbability[n] = I.Fitness;
        estimatedProbabilitySum += std::exp(I.Fitness);

        if (I.Fitness > maxValue)
            maxValue = I.Fitness;
    }


    Eigen::VectorXd normalisedProbability = (estimatedProbability - maxValue * Eigen::VectorXd::Ones(N)).array().exp();
    normalisedProbability = normalisedProbability / normalisedProbability.sum();
    return Rcpp::List::create(Rcpp::Named("GenotypeMatrix") = genotypeMatrices,
                              Rcpp::Named("LogUnnormalisedProbabilitySum") = std::log(estimatedProbabilitySum),
                              Rcpp::Named("LogUnnormalisedProbability") = estimatedProbability,
                              Rcpp::Named("NormalisedProbabilities") = normalisedProbability);
}


////
struct sampledPriorGenotype {
    Eigen::VectorXd encodedGenotype;
    double priorProbability;

    sampledPriorGenotype(Eigen::VectorXd encodedGenotype_, double priorProbability_) :
        encodedGenotype(encodedGenotype_), priorProbability(priorProbability_) {};
};


double calculatePriorProbabiliy(const Eigen::VectorXd & E, const Eigen::MatrixXd &knownProfiles, const Eigen::VectorXd &coverage,
                                const double &sampleParameter, const Eigen::VectorXd &mixtureParameters, const double &markerParameter,
                                const std::size_t & numberOfContributors, const bool &suggestionBool)
{
    const std::size_t numberOfKnownContributors = knownProfiles.cols();
    const std::size_t N = E.size();
    const std::size_t M = coverage.size();

    Eigen::MatrixXd genotype = Eigen::MatrixXd::Zero(M, numberOfContributors);
    for (std::size_t k = 0; k < numberOfKnownContributors; k++)
    {
        genotype.col(k) = knownProfiles.col(k);
    }

    double priorProbability = 0.0;
    Eigen::VectorXd adjustedCoverage = Eigen::VectorXd::Ones(M);
    for (std::size_t n = 0; n < N; n++)
    {
        if (suggestionBool)
        {
            adjustedCoverage = coverage - (markerParameter * sampleParameter) * (genotype * mixtureParameters);
            for (std::size_t m = 0; m < M; m++)
            {
                if (adjustedCoverage[m] < 0.0)
                {
                    adjustedCoverage[m] = 1.0;
                }
            }
        }

        std::size_t u = std::floor(n / 2);
        Eigen::VectorXd proportionalProbability = adjustedCoverage / adjustedCoverage.sum();

        genotype(E[n], numberOfKnownContributors + u) += 1.0;
        priorProbability += std::log(proportionalProbability[E[n]]);

    }

    return priorProbability;
}

sampledPriorGenotype samplePriorGenotypeMarker(const Eigen::MatrixXd &knownProfiles, const Eigen::VectorXd &coverage,
                                               const double &sampleParameter, const Eigen::VectorXd &mixtureParameters,
                                               const double &markerParameter, const std::size_t numberOfContributors, const bool & suggestionBool,
                                               const std::size_t &seed)
{
    boost::random::mt19937 rng(seed);
    boost::random::uniform_real_distribution<> uniform(0, 1);

    const std::size_t &N = coverage.size();
    const std::size_t &numberOfKnownContributors = knownProfiles.cols();
    const std::size_t &numberOfUnknownContributors = numberOfContributors - numberOfKnownContributors;
    Eigen::VectorXd adjustedCoverage = Eigen::VectorXd::Ones(N);

    Eigen::MatrixXd sampledGenotype = Eigen::MatrixXd::Zero(N, numberOfContributors);
    for (std::size_t k = 0; k < numberOfKnownContributors; k++)
        sampledGenotype.col(k) = knownProfiles.col(k);

    Eigen::VectorXd sampledGenotypeEncoded = Eigen::VectorXd::Zero(2 * numberOfUnknownContributors);
    double sampledPriorProbability = 0.0;
    for (std::size_t u = 0; u < numberOfUnknownContributors; u++)
    {
        for (std::size_t i = 0; i < 2; i++)
        {
            if (suggestionBool)
            {
                adjustedCoverage = coverage - (markerParameter * sampleParameter) * (sampledGenotype * mixtureParameters);
                for (std::size_t n = 0; n < N; n++)
                {
                    if (adjustedCoverage[n] < 0.0)
                    {
                        adjustedCoverage[n] = 1.0;
                    }
                }
            }

            Eigen::VectorXd proportionalProbability = adjustedCoverage / adjustedCoverage.sum();
            Eigen::VectorXd partialSumProportionalProbability = partialSumEigen(proportionalProbability);
            double randomVariate = uniform(rng);
            std::size_t j = 0.0;
            while (randomVariate > partialSumProportionalProbability[j + 1])
                j++;

            sampledGenotype(j, numberOfKnownContributors + u) += 1.0;
            sampledPriorProbability += std::log(proportionalProbability[j]);

            if (i == 0)
            {
                sampledGenotypeEncoded[2 * u] = j;
            }
            else
            {
                if (j < sampledGenotypeEncoded[2 * u])
                {
                    sampledGenotypeEncoded[2 * u + i] = sampledGenotypeEncoded[2 * u];
                    sampledGenotypeEncoded[2 * u] = j;
                }
                else
                {
                    sampledGenotypeEncoded[2 * u + i] = j;
                }
            }
        }
    }

    sampledPriorGenotype spg(sampledGenotypeEncoded, sampledPriorProbability);
    return spg;
}

//[[Rcpp::export(.samplePosteriorGenotypesGuidedCpp)]]
std::vector<Eigen::MatrixXd> samplePosteriorGenotypesGuidedCpp(const Eigen::VectorXd & encodedProfiles, const Eigen::VectorXd & sampleParameters,
                                                               const Eigen::VectorXd & noiseParameters, const Eigen::VectorXd & mixtureParameters,
                                                               const Eigen::VectorXd markerParameters, const Eigen::VectorXd & coverage,
                                                               const std::vector< std::vector< Eigen::MatrixXd > > potentialParents,
                                                               const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                                                               const Eigen::VectorXd & alleleFrequencies, const double & theta,
                                                               const std::size_t & numberOfContributors,  const Eigen::VectorXd & numberOfAlleles,
                                                               const std::size_t & levelsOfStutterRecursion,
                                                               const std::size_t &numberOfSimulations, const bool & suggestionBool,
                                                               const std::size_t &seed)
{
    boost::random::mt19937 rng(seed);
    boost::random::uniform_int_distribution<> seedShift(0, seed);
    boost::random::uniform_real_distribution<> uniform(0, 1);

    const std::size_t &numberOfKnownContributors = knownProfiles.cols();
    const std::size_t &numberOfMarkers = numberOfAlleles.size();
    const Eigen::VectorXd &partialSumAlleles = partialSumEigen(numberOfAlleles);

    const ExperimentalSetup ES(1, numberOfAlleles, numberOfContributors, numberOfKnownContributors,
                               knownProfiles, allKnownProfiles, coverage, potentialParents, markerParameters, 0.0,
                               0.0, theta, alleleFrequencies, levelsOfStutterRecursion);

    Eigen::VectorXd currentEncodedGenotype = encodedProfiles;
    double currentPriorProbability = calculatePriorProbabiliy(currentEncodedGenotype, knownProfiles, coverage, sampleParameters[0], mixtureParameters,
                                                              markerParameters[0], numberOfContributors, suggestionBool);

    Individual currentI(currentEncodedGenotype, sampleParameters, noiseParameters, mixtureParameters, markerParameters, ES);
    double currentLogLikelihood = currentI.Fitness;

    std::vector<Eigen::MatrixXd> sampledGenotypeList(numberOfSimulations);
    for (std::size_t n = 0; n < numberOfSimulations; n++)
    {
        int seedShift_n = seedShift(rng);
        sampledPriorGenotype newGenotype = samplePriorGenotypeMarker(knownProfiles, coverage, sampleParameters[0], mixtureParameters,
                                                                     markerParameters[0], numberOfContributors, suggestionBool,
                                                                     seed + n + seedShift_n);

        Eigen::VectorXd newEncodedGenotype = newGenotype.encodedGenotype;
        double newPriorProbability = newGenotype.priorProbability;

        Individual newI(newEncodedGenotype, sampleParameters, noiseParameters, mixtureParameters, markerParameters, ES);
        double newLogLikelihood = newI.Fitness;

        double logHastingsRatio = newLogLikelihood - currentLogLikelihood + currentPriorProbability - newPriorProbability;
        double acceptProbability = 1.0;
        if (logHastingsRatio < 0)
        {
            acceptProbability = std::exp(logHastingsRatio);
        }

        double randomVariate = uniform(rng);
        if (randomVariate < acceptProbability)
        {
            currentEncodedGenotype = newEncodedGenotype;
            currentPriorProbability = newPriorProbability;
            currentLogLikelihood = newLogLikelihood;
        }

        Eigen::MatrixXd currentDecodedGenotype = decoding(currentEncodedGenotype, numberOfAlleles, 1, numberOfContributors - numberOfKnownContributors);
        sampledGenotypeList[n] = currentDecodedGenotype;
    }

    return sampledGenotypeList;
}

