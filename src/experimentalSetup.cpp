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


InitialPopulationRandomVariates::InitialPopulationRandomVariates(const Eigen::VectorXd & numberOfAlleles, const std::size_t & seed) :
    rng(seed), M(numberOfAlleles.size())
{
    for (std::size_t m = 0; m < M; m++)
    {
        boost::random::uniform_int_distribution<> uniform_interval_m(0, numberOfAlleles[m] - 1);
        interval_uniform_generator uniform_interval_vector_m(rng, uniform_interval_m);
        generate_uniform_interval.push_back(uniform_interval_vector_m);
    }
}


RandomVariates::RandomVariates(const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfUnknownContributors,
                               const std::size_t & seed) :
    rng(seed), M(numberOfAlleles.size()),
    uniform_real(0.0, 1.0), generate_uniform_real(rng, uniform_real),
    uniform_01(0, 1), generate_uniform_binary(rng, uniform_01),
    uniform_marker(0, M - 1), generate_uniform_marker(rng, uniform_marker),
    uniform_unknown_contributor(0, numberOfUnknownContributors - 1), generate_uniform_unknown_contributor(rng, uniform_unknown_contributor)
{
    for (std::size_t m = 0; m < M; m++)
    {
        boost::random::uniform_int_distribution<> uniform_mutation_m(1, numberOfAlleles[m] - 1);
        interval_uniform_generator uniform_mutation_vector_m(rng, uniform_mutation_m);
        generate_uniform_mutation.push_back(uniform_mutation_vector_m);
    }
}

//[[Rcpp::export(.testing_rv_struct)]]
Eigen::MatrixXd testing_rv_struct(const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfUnknownContributors,
                                  const std::size_t & seed, const std::size_t & numberOfSimulations)
{
    RandomVariates RV(numberOfAlleles, numberOfUnknownContributors, seed);

    const std::size_t & N = numberOfSimulations;
    const std::size_t & M = numberOfAlleles.size();

    Eigen::MatrixXd res(N, M + 4);
    for (std::size_t n = 0; n < N; n++)
    {
        res(n, 0) = RV.generate_uniform_real();
        res(n, 1) = RV.generate_uniform_binary();
        res(n, 2) = RV.generate_uniform_marker();
        res(n, 3) = RV.generate_uniform_unknown_contributor();
        for (std::size_t m = 0; m < M; m++)
        {
            res(n, 4 + m) = RV.generate_uniform_mutation[m]();
        }
    }

    return res;
}

//
ExperimentalSetup::ExperimentalSetup(const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfContributors,
                                     const std::size_t & numberOfKnownContributors, const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                                     const Eigen::VectorXd & coverage, const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents, const Eigen::VectorXd & markerImbalances,
                                     const double & convexMarkerImbalanceInterpolation,
                                     const Eigen::VectorXd & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies, const std::size_t & levelsOfStutterRecursion)
{
    NumberOfMarkers = numberOfMarkers;
    NumberOfAlleles = numberOfAlleles;
    PartialSumAlleles = partialSumEigen(numberOfAlleles);

    NumberOfContributors = numberOfContributors;
    KnownProfiles = knownProfiles;
    NumberOfKnownContributors = numberOfKnownContributors;

    AllKnownProfiles = allKnownProfiles;

    Coverage = coverage;
    PotentialParents = potentialParents;
    MarkerImbalances = markerImbalances;

    ConvexMarkerImbalanceInterpolation = convexMarkerImbalanceInterpolation;
    Tolerance = tolerance;
    LevelsOfStutterRecursion = levelsOfStutterRecursion;

    Theta = theta;
    AlleleFrequencies = alleleFrequencies;
}


Eigen::VectorXd ExperimentalSetup::GenerateUnknownGenotype(InitialPopulationRandomVariates & RV)
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

            double d = RV.generate_uniform_interval[i]();

            individual(k) = d;
        }
    }

    return individual;
}

