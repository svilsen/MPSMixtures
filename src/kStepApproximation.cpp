#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include "experimentalSetup.hpp"
#include "individual.hpp"
#include "encodingScheme.hpp"
#include "logLikelihoods.hpp"

#include "AuxiliaryFunctions.hpp"

//[[Rcpp::export(.kStepApproximation)]]
double kStepApproximation(const Eigen::VectorXd & encodedProfiles, const Eigen::VectorXd & sampleParameters,
                          const Eigen::VectorXd & noiseParameters, const Eigen::VectorXd & mixtureParameters,
                          const Eigen::VectorXd & coverage, const Eigen::VectorXd markerImbalances,
                          const std::vector< std::vector< Eigen::MatrixXd > > potentialParents,
                          const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                          const Eigen::VectorXd & alleleFrequencies, const double & theta,
                          const std::size_t & numberOfContributors, const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles,
                          const std::size_t & levelsOfStutterRecursion,
                          const std::size_t & kStep, const std::size_t & markerIndicator, const double & normalisingConstant)
{
    const std::size_t numberOfKnownContributors = knownProfiles.cols();
    const ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors,
                               knownProfiles, allKnownProfiles, coverage, potentialParents, markerImbalances, 0.8,
                               1e-6, theta, alleleFrequencies, levelsOfStutterRecursion);

    const Individual optimal(encodedProfiles, sampleParameters, noiseParameters, mixtureParameters, markerImbalances, ES);
    const std::size_t & numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;

    double estimatedProbability = 0.0;
    for (std::size_t u = 0; u < numberOfUnknownContributors; u++)
    {
        for (std::size_t m = markerIndicator; m < numberOfMarkers; m++)
        {
            Eigen::VectorXd encodedProfiles_um = encodedProfiles;
            std::size_t a = 2 * m * numberOfUnknownContributors + 2 * u;
            for (std::size_t i = 0; i < (numberOfAlleles[m] - 1); i++)
            {
                Eigen::VectorXd encodedProfiles_um = encodedProfiles;
                encodedProfiles_um[a] = static_cast<int>(encodedProfiles_um[a] + 1) % static_cast<int>(numberOfAlleles[m]);

                for (std::size_t j = 0; j < (numberOfAlleles[m] - 1); j++)
                {
                    encodedProfiles_um[a + 1] = (encodedProfiles_um[a + 1] + 1);
                    Individual I_umij(encodedProfiles_um, sampleParameters, noiseParameters, mixtureParameters, markerImbalances, ES);
                    estimatedProbability += std::exp(I_umij.Fitness - normalisingConstant);

                    if (kStep > 0)
                    {
                        estimatedProbability += kStepApproximation(encodedProfiles, sampleParameters, noiseParameters, mixtureParameters,
                                                                   coverage, markerImbalances, potentialParents,
                                                                   knownProfiles, allKnownProfiles, alleleFrequencies, theta,
                                                                   numberOfContributors, numberOfMarkers, numberOfAlleles, levelsOfStutterRecursion,
                                                                   kStep - 1, m + 1, normalisingConstant);
                    }
                }
            }
        }
    }

    return estimatedProbability;
}
