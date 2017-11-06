#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//[[Rcpp::depends(BH)]]
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "experimentalSetup.hpp"
#include "encodingScheme.hpp"
#include "AuxiliaryFunctions.hpp"

ExperimentalSetup::ExperimentalSetup(std::size_t numberOfMarkers, Eigen::VectorXd numberOfAlleles, std::size_t numberOfContributors, std::size_t numberOfKnownContributors,
                                     Eigen::MatrixXd knownProfiles, Eigen::MatrixXd allKnownProfiles,
                                     Eigen::VectorXd coverage, std::vector< std::vector < Eigen::MatrixXd > > potentialParents, Eigen::VectorXd markerImbalances,
                                     double tolerance, double theta, Eigen::VectorXd alleleFrequencies)
{
    NumberOfMarkers = numberOfMarkers;
    NumberOfAlleles = numberOfAlleles;

    NumberOfContributors = numberOfContributors;
    KnownProfiles = knownProfiles;
    NumberOfKnownContributors = numberOfKnownContributors;

    AllKnownProfiles = allKnownProfiles;

    Coverage = coverage;
    PotentialParents = potentialParents;
    MarkerImbalances = markerImbalances;

    Tolerance = tolerance;
    Theta = theta;
    AlleleFrequencies = alleleFrequencies;
}

Eigen::VectorXd ExperimentalSetup::GenerateUnknownGenotype(std::size_t seed)
{
    const Eigen::VectorXd partialSumAlleles = partialSumEigen(NumberOfAlleles);
    std::size_t NumberOfUnknownContributors = (NumberOfContributors - NumberOfKnownContributors);

    Eigen::VectorXd individual(2 * NumberOfUnknownContributors * NumberOfMarkers);
    for (std::size_t i = 0; i < NumberOfMarkers; i++)
    {

        boost::random::uniform_int_distribution<> uniform(0, NumberOfAlleles[i] - 1);
        for (std::size_t j = 0; j < 2 * NumberOfUnknownContributors; j++)
        {
            std::size_t k = 2 * NumberOfUnknownContributors * i + j;
            boost::random::mt19937 rng(seed + k);
            double d = uniform(rng);

            individual(k) = d;
        }
    }

    return individual;
}

