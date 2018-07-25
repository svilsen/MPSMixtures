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
                                                        const double & convexMarkerImbalanceInterpolation, const Eigen::VectorXd & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                                                        const std::size_t & populationSize, const std::size_t & numberOfIterations, const std::size_t & numberOfIterationsEqualMinMax,
                                                        const std::size_t & numberOfFittestIndividuals, const int & parentSelectionWindowSize, const bool & allowParentSurvival,
                                                        const double & fractionEnsuredSurvival, const double & crossoverProbability, const double & mutationProbabilityLowerLimit,
                                                        const std::size_t & mutationIterations, const double & mutationDegreesOfFreedom, const Eigen::VectorXd mutationDecay,
                                                        const std::size_t hillClimbingIterations,
                                                        const std::size_t & seed, const bool & trace, const std::size_t & levelsOfStutterRecursion)
{
    ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                         coverage, potentialParents, markerImbalances, convexMarkerImbalanceInterpolation, tolerance, theta, alleleFrequencies, levelsOfStutterRecursion);

    InitialPopulationRandomVariates InitialRV(numberOfAlleles, seed);

    EvolutionaryAlgorithm EA(ES, populationSize, numberOfIterations, numberOfIterationsEqualMinMax, numberOfFittestIndividuals, parentSelectionWindowSize, allowParentSurvival,
                             crossoverProbability, mutationProbabilityLowerLimit, mutationIterations, mutationDegreesOfFreedom, mutationDecay, fractionEnsuredSurvival,
                             hillClimbingIterations, InitialRV);

    RandomVariates RV(numberOfAlleles, numberOfContributors - numberOfKnownContributors, seed + 1);
    EA.Run(ES, RV, trace);

    Population FittestIndividuals = EA.FittestMembersOfEntireRun;
    Rcpp::List RL = FittestIndividuals.ReturnRcppList(ES);

    return Rcpp::List::create(Rcpp::Named("U") = RL,
                              Rcpp::Named("UMaxSize") = numberOfFittestIndividuals,
                              Rcpp::Named("NumberOfIterations") = EA.Iteration,
                              Rcpp::Named("NumberOfIterationsEqualMax") = numberOfIterationsEqualMinMax,
                              Rcpp::Named("NumberOfPopulations") = 1);
}

//[[Rcpp::export(.initialisingParallelEvolutionaryAlgorithm)]]
Rcpp::List initialisingParallelEvolutionaryAlgorithm(const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfContributors,
                                                     const std::size_t & numberOfKnownContributors, const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                                                     const Eigen::VectorXd & coverage, const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents, const Eigen::VectorXd & markerImbalances,
                                                     const double & convexMarkerImbalanceInterpolation, const Eigen::VectorXd & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                                                     const std::size_t & populationSize, const std::size_t & seed, const std::size_t & levelsOfStutterRecursion)
{
    ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                         coverage, potentialParents, markerImbalances, convexMarkerImbalanceInterpolation, tolerance, theta, alleleFrequencies, levelsOfStutterRecursion);

    InitialPopulationRandomVariates InitialRV(numberOfAlleles, seed);
    EvolutionaryAlgorithm EA(ES, populationSize, InitialRV);

    Rcpp::List RL = EA.CurrentPopulation.ReturnCompressedRcppList();
    return RL;
}

//[[Rcpp::export(.runningParallelEvolutionaryAlgorithm)]]
Rcpp::List runningParallelEvolutionaryAlgorithm(const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfContributors, const std::size_t & numberOfKnownContributors,
                                                const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                                                const Eigen::VectorXd & coverage, const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents, const Eigen::VectorXd & markerImbalances,
                                                const double & convexMarkerImbalanceInterpolation,
                                                const Eigen::VectorXd & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                                                const std::size_t & numberOfIterations, const std::size_t & numberOfIterationsEqualMinMax,
                                                const std::size_t & numberOfFittestIndividuals, const int & parentSelectionWindowSize, const bool allowParentSurvival,
                                                const double & fractionEnsuredSurvival, const double & crossoverProbability, const double & mutationProbabilityLowerLimit,
                                                const std::size_t & mutationIterations, const double & mutationDegreesOfFreedom, const Eigen::VectorXd & mutationDecay,
                                                const std::size_t & hillClimbingIterations, const std::size_t & seed, const bool & trace,
                                                const Eigen::MatrixXd encodedPopulationList, const Eigen::MatrixXd sampleParametersList,
                                                const Eigen::MatrixXd noiseParametersList, const Eigen::MatrixXd mixtureParametersList,
                                                const Eigen::MatrixXd markerParametersList,
                                                const Eigen::VectorXd fitnessList, const std::size_t & levelsOfStutterRecursion)
{
    ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                         coverage, potentialParents, markerImbalances, convexMarkerImbalanceInterpolation, tolerance, theta, alleleFrequencies, levelsOfStutterRecursion);

    Population currentPopulation(encodedPopulationList, sampleParametersList, noiseParametersList, mixtureParametersList, markerParametersList, fitnessList, ES);

    EvolutionaryAlgorithm EA(ES, currentPopulation,
                             numberOfIterations, numberOfIterationsEqualMinMax, numberOfFittestIndividuals,
                             parentSelectionWindowSize, allowParentSurvival,
                             crossoverProbability, mutationProbabilityLowerLimit, mutationIterations,
                             mutationDegreesOfFreedom, mutationDecay, fractionEnsuredSurvival, hillClimbingIterations);


    RandomVariates RV(numberOfAlleles, numberOfContributors - numberOfKnownContributors, seed);
    EA.Run(ES, RV, trace);

    Population FittestIndividuals = EA.FittestMembersOfEntireRun;

    Rcpp::List RL = EA.CurrentPopulation.ReturnCompressedRcppList();
    Rcpp::List FL = FittestIndividuals.ReturnRcppList(ES);
    return Rcpp::List::create(Rcpp::Named("LastPopulation") = RL,
                              Rcpp::Named("FittestIndividuals") = FL);
}
