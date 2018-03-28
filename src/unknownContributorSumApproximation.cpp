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

    const Eigen::VectorXd tolerance = Eigen::VectorXd::Zero(4);
    const ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors,
                               knownProfiles, allKnownProfiles, coverage, potentialParents, markerImbalances, 0.8,
                               tolerance, theta, alleleFrequencies, levelsOfStutterRecursion);

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

    const Eigen::VectorXd tolerance = Eigen::VectorXd::Zero(4);
    const ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors,
                               knownProfiles, allKnownProfiles, coverage, potentialParents, markerImbalances, 0.8,
                               tolerance, theta, alleleFrequencies, levelsOfStutterRecursion);

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

Eigen::MatrixXd expectedContributionMarker(const Eigen::MatrixXd & genotype, const std::vector<Eigen::MatrixXd> &potentialParents_m,
                                           const double & levelsOfStutterRecursion)
{
    const std::size_t & N = genotype.rows();
    const std::size_t & M = genotype.cols();

    Eigen::MatrixXd expectedContributionProfile_m = Eigen::MatrixXd::Zero(N, M);
    for (std::size_t u = 0; u < M; u++)
    {
        Eigen::VectorXd genotype_u = genotype.col(u);
        Eigen::VectorXd stutterContribution = Eigen::VectorXd::Zero(N);
        for (std::size_t n = 0; n < N; n++)
        {
            stutterContribution[n] = ParentStutterContribution(n, levelsOfStutterRecursion, 1, genotype_u, potentialParents_m, N);
        }

        expectedContributionProfile_m.col(u) = genotype_u + stutterContribution;
    }

    return expectedContributionProfile_m;
}


void adjustCoverage(Eigen::VectorXd & adjustedCoverage, const double & sampleParameter, const Eigen::VectorXd & mixtureParamters,
                    const double & markerParameter, const Eigen::MatrixXd & genotype, const std::vector< Eigen::MatrixXd > & potentialParents_m,
                    const double & levelsOfStutterRecursion)
{
    const std::size_t & N = adjustedCoverage.size();
    const Eigen::VectorXd & E_c = expectedContributionMarker(genotype, potentialParents_m, levelsOfStutterRecursion);
    const Eigen::VectorXd & mu_ma = sampleParameter * markerParameter * (E_c * mixtureParamters);

    adjustedCoverage = adjustedCoverage - mu_ma;
    for (std::size_t n = 0; n < N; n++)
    {
        if (adjustedCoverage[n] <= 0.0)
        {
            adjustedCoverage[n] = std::exp(adjustedCoverage[n]);
        }

    }
}

double calculatePriorProbabiliy(const Eigen::VectorXd & E, const Eigen::MatrixXd &knownProfiles, const Eigen::VectorXd &coverage,
                                std::vector<Eigen::MatrixXd> potentialParents_m,
                                const double &sampleParameter, const Eigen::VectorXd &mixtureParameters, const double &markerParameter,
                                const std::size_t & numberOfContributors, const double & levelsOfStutterRecursion,
                                const bool &suggestionBool)
{
    if (!suggestionBool)
    {
        return 0.0;
    }

    const std::size_t numberOfKnownContributors = knownProfiles.cols();
    const std::size_t N = E.size();
    const std::size_t M = coverage.size();

    Eigen::MatrixXd genotype = Eigen::MatrixXd::Zero(M, numberOfContributors);
    for (std::size_t k = 0; k < numberOfKnownContributors; k++)
    {
        genotype.col(k) = knownProfiles.col(k);
    }

    Eigen::VectorXd adjustedCoverage = coverage;
    if (suggestionBool)
    {
        adjustCoverage(adjustedCoverage, sampleParameter, mixtureParameters, markerParameter, genotype, potentialParents_m, levelsOfStutterRecursion);
    }
    else
    {
        adjustedCoverage = Eigen::VectorXd::Ones(M);
    }

    double priorProbability = 0.0;
    for (std::size_t n = 0; n < N; n++)
    {
        std::size_t u = std::floor(n / 2);
        Eigen::VectorXd proportionalProbability = adjustedCoverage / adjustedCoverage.sum();

        genotype(E[n], numberOfKnownContributors + u) += 1.0;

        if (suggestionBool)
        {
            Eigen::VectorXd singleAllele = Eigen::VectorXd::Zero(M);
            singleAllele[E[n]] = 1.0;

            priorProbability += std::log(proportionalProbability[E[n]]);
            adjustCoverage(adjustedCoverage, sampleParameter, mixtureParameters.row(u), markerParameter, singleAllele, potentialParents_m,
                           levelsOfStutterRecursion);
        }
    }

    return priorProbability;
}

sampledPriorGenotype samplePriorGenotypeMarker(const Eigen::MatrixXd &knownProfiles, const Eigen::VectorXd &coverage,
                                               const std::vector<Eigen::MatrixXd> &potentialParents_m,
                                               const double &sampleParameter, const Eigen::VectorXd &mixtureParameters,
                                               const double &markerParameter, const std::size_t numberOfContributors,
                                               const double &levelsOfStutterRecursion, const bool & suggestionBool,
                                               const std::size_t &seed)
{
    boost::random::mt19937 rng(seed);
    boost::random::uniform_real_distribution<> uniform(0, 1);

    const std::size_t &N = coverage.size();
    const std::size_t &numberOfKnownContributors = knownProfiles.cols();
    const std::size_t &numberOfUnknownContributors = numberOfContributors - numberOfKnownContributors;

    Eigen::MatrixXd sampledGenotype = Eigen::MatrixXd::Zero(N, numberOfContributors);
    for (std::size_t k = 0; k < numberOfKnownContributors; k++)
        sampledGenotype.col(k) = knownProfiles.col(k);

    Eigen::VectorXd adjustedCoverage = Eigen::VectorXd::Ones(N);
    if (suggestionBool)
    {
        adjustedCoverage = coverage;
        adjustCoverage(adjustedCoverage, sampleParameter, mixtureParameters, markerParameter, sampledGenotype, potentialParents_m,
                       levelsOfStutterRecursion);
    }

    Eigen::VectorXd sampledGenotypeEncoded = Eigen::VectorXd::Zero(2 * numberOfUnknownContributors);
    double sampledPriorProbability = 0.0;
    for (std::size_t u = 0; u < numberOfUnknownContributors; u++)
    {
        for (std::size_t i = 0; i < 2; i++)
        {
            Eigen::VectorXd proportionalProbability = adjustedCoverage / adjustedCoverage.sum();
            Eigen::VectorXd partialSumProportionalProbability = partialSumEigen(proportionalProbability);
            double randomVariate = uniform(rng);
            std::size_t j = 0.0;
            while (randomVariate > partialSumProportionalProbability[j + 1])
                j++;

            sampledGenotype(j, numberOfKnownContributors + u) += 1.0;

            if (suggestionBool)
            {
                Eigen::VectorXd singleAllele = Eigen::VectorXd::Zero(N);
                singleAllele[j] = 1.0;

                sampledPriorProbability += std::log(proportionalProbability[j]);
                adjustCoverage(adjustedCoverage, sampleParameter, mixtureParameters.row(u), markerParameter, singleAllele, potentialParents_m,
                               levelsOfStutterRecursion);
            }

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
Rcpp::List samplePosteriorGenotypesGuidedCpp(const Eigen::VectorXd & encodedProfiles, const Eigen::VectorXd & sampleParameters,
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

    const Eigen::VectorXd tolerance = Eigen::VectorXd::Zero(4);
    const ExperimentalSetup ES(1, numberOfAlleles, numberOfContributors, numberOfKnownContributors,
                               knownProfiles, allKnownProfiles, coverage, potentialParents, markerParameters, 0.0,
                               tolerance, theta, alleleFrequencies, levelsOfStutterRecursion);

    Eigen::VectorXd currentEncodedGenotype = encodedProfiles;
    double currentPriorProbability = calculatePriorProbabiliy(currentEncodedGenotype, knownProfiles, coverage,
                                                              potentialParents[0], sampleParameters[0], mixtureParameters,
                                                              markerParameters[0], numberOfContributors, levelsOfStutterRecursion,
                                                              suggestionBool);

    Individual currentI(currentEncodedGenotype, sampleParameters, noiseParameters, mixtureParameters, markerParameters, ES);
    double currentLogLikelihood = currentI.Fitness;

    std::vector<Eigen::MatrixXd> sampledGenotypeList(numberOfSimulations);
    double acceptedProposals = 0.0;
    for (std::size_t n = 0; n < numberOfSimulations; n++)
    {
        int seedShift_n = seedShift(rng);
        const sampledPriorGenotype & newGenotype = samplePriorGenotypeMarker(knownProfiles, coverage, potentialParents[0],
                                                                             sampleParameters[0], mixtureParameters,
                                                                             markerParameters[0], numberOfContributors, levelsOfStutterRecursion,
                                                                             suggestionBool, seed + n + seedShift_n);

        const Eigen::VectorXd & newEncodedGenotype = newGenotype.encodedGenotype;
        const double & newPriorProbability = newGenotype.priorProbability;

        const Individual newI(newEncodedGenotype, sampleParameters, noiseParameters, mixtureParameters, markerParameters, ES);
        const double & newLogLikelihood = newI.Fitness;

        const double & logHastingsRatio = newLogLikelihood - currentLogLikelihood + currentPriorProbability - newPriorProbability;
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
            acceptedProposals += 1;
        }

        Eigen::MatrixXd currentDecodedGenotype = decoding(currentEncodedGenotype, numberOfAlleles, 1,
                                                          numberOfContributors - numberOfKnownContributors);
        sampledGenotypeList[n] = currentDecodedGenotype;
    }

    return Rcpp::List::create(Rcpp::Named("SampledGenotypes") = sampledGenotypeList,
                              Rcpp::Named("AcceptedProposals") = acceptedProposals);
}
