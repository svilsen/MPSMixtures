#include <Rcpp.h>

#include <RcppEigen.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "evolutionaryAlgorithm.hpp"
#include "AuxiliaryFunctions.hpp"

//[[Rcpp::export(.runningSinglePopulationEvolutionaryAlgorithm)]]
Rcpp::List runningSinglePopulationEvolutionaryAlgorithm(const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfContributors, const std::size_t & numberOfKnownContributors,
                                                        const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                                                        const Eigen::VectorXd & coverage, const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents, const Eigen::VectorXd & markerImbalances,
                                                        const double & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                                                        const std::size_t & populationSize, const std::size_t & numberOfIterations, const std::size_t & numberOfIterationsEqualMinMax,
                                                        const std::size_t & numberOfFittestIndividuals, const int & parentSelectionWindowSize, const bool & allowParentSurvival,
                                                        const double & crossoverProbability, const double & mutationProbabilityLowerLimit, const double & mutationDegreesOfFreedom,
                                                        const Eigen::VectorXd mutationDecay,
                                                        const std::size_t hillClimbingDirections, const std::size_t hillClimbingIterations,
                                                        const std::size_t & seed, const bool & trace)
{
    ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                         coverage, potentialParents, markerImbalances, tolerance, theta, alleleFrequencies);

    boost::random::mt19937 rngSeed(seed);
    boost::random::uniform_int_distribution<> uniformShift(0, 1e6);

    std::size_t seedShift = uniformShift(rngSeed);
    EvolutionaryAlgorithm EA(ES, populationSize, numberOfIterations, numberOfIterationsEqualMinMax, numberOfFittestIndividuals, parentSelectionWindowSize, allowParentSurvival,
                             crossoverProbability, mutationProbabilityLowerLimit, mutationDegreesOfFreedom, mutationDecay, 1,
                             hillClimbingDirections, hillClimbingIterations, seed + seedShift);

    std::size_t seedShift2 = uniformShift(rngSeed);

    EA.Run(ES, seed + seedShift, trace);

    Population FittestIndividuals = EA.FittestMembersOfEntireRun;

    Rcpp::List RL = FittestIndividuals.ReturnRcppList();
    return RL;
}


//[[Rcpp::export(.initialisingParallelEvolutionaryAlgorithm)]]
Rcpp::List initialisingParallelEvolutionaryAlgorithm(const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfContributors,
                                                     const std::size_t & numberOfKnownContributors, const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                                                     const Eigen::VectorXd & coverage, const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents, const Eigen::VectorXd & markerImbalances,
                                                     const double & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                                                     const std::size_t populationSize, const std::size_t seed)
{
    ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                         coverage, potentialParents, markerImbalances, tolerance, theta, alleleFrequencies);

    boost::random::mt19937 rngSeed(seed);
    boost::random::uniform_int_distribution<> uniformShift(0, 1e6);

    std::size_t seedShift = uniformShift(rngSeed);
    EvolutionaryAlgorithm EA(ES, populationSize, seed + seedShift);

    Rcpp::List RL = EA.CurrentPopulation.ReturnCompressedRcppList();
    return RL;
}

//[[Rcpp::export(.runningParallelEvolutionaryAlgorithm)]]
Rcpp::List runningParallelEvolutionaryAlgorithm(const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfContributors, const std::size_t & numberOfKnownContributors,
                                                const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                                                const Eigen::VectorXd & coverage, const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents, const Eigen::VectorXd & markerImbalances,
                                                const double & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                                                const std::size_t & numberOfIterations, const std::size_t & numberOfIterationsEqualMinMax,
                                                const std::size_t & numberOfFittestIndividuals, const int & parentSelectionWindowSize, const bool & allowParentSurvival,
                                                const double & crossoverProbability, const double & mutationProbabilityLowerLimit, const double & mutationDegreesOfFreedom,
                                                const Eigen::VectorXd mutationDecay,
                                                const std::size_t hillClimbingDirections, const std::size_t hillClimbingIterations,
                                                const std::size_t & seed, const bool & trace,
                                                const Eigen::MatrixXd encodedPopulationList, const Eigen::MatrixXd sampleParametersList,
                                                const Eigen::MatrixXd noiseParametersList, const Eigen::MatrixXd mixtureParametersList)
{
    ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                         coverage, potentialParents, markerImbalances, tolerance, theta, alleleFrequencies);
    const std::size_t populationSize = encodedPopulationList.cols();

    std::vector< Individual > currentIndividuals(populationSize);
    for (std::size_t i = 0; i < populationSize; i++)
    {
        Individual I(encodedPopulationList.col(i), sampleParametersList.col(i), noiseParametersList.col(i), mixtureParametersList.col(i), ES);
        currentIndividuals[i] = I;
    }

    Population currentPopulation(currentIndividuals);

    EvolutionaryAlgorithm EA(ES, currentPopulation,
                             numberOfIterations, numberOfIterationsEqualMinMax, numberOfFittestIndividuals,
                             parentSelectionWindowSize, allowParentSurvival,
                             crossoverProbability, mutationProbabilityLowerLimit, mutationDegreesOfFreedom, mutationDecay, 1,
                             hillClimbingDirections, hillClimbingIterations);

    boost::random::mt19937 rngSeed(seed);
    boost::random::uniform_int_distribution<> uniformShift(0, 1e6);

    std::size_t seedShift = uniformShift(rngSeed);

    EA.Run(ES, seed + seedShift, trace);

    Rcpp::List RL = EA.CurrentPopulation.ReturnCompressedRcppList();
    Rcpp::List FL = EA.FittestMembersOfEntireRun.ReturnRcppList();
    return Rcpp::List::create(Rcpp::Named("LastPopulation") = RL,
                              Rcpp::Named("FittestIndividuals") = FL);
}
