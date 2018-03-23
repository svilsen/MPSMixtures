// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// partialSumEigen
Eigen::VectorXd partialSumEigen(const Eigen::VectorXd& x);
RcppExport SEXP _MPSMixtures_partialSumEigen(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(partialSumEigen(x));
    return rcpp_result_gen;
END_RCPP
}
// setupIndividual
Rcpp::List setupIndividual(const std::size_t& numberOfMarkers, const Eigen::VectorXd& numberOfAlleles, const std::size_t& numberOfContributors, const std::size_t& numberOfKnownContributors, const Eigen::MatrixXd& knownProfiles, const Eigen::VectorXd& coverage, const std::vector< std::vector < Eigen::MatrixXd > >& potentialParents, const Eigen::VectorXd& markerImbalances, const double& convexMarkerImbalanceInterpolation, const double& tolerance, const double& theta, const Eigen::VectorXd& alleleFrequencies, const std::size_t& levelsOfStutterRecursion);
RcppExport SEXP _MPSMixtures_setupIndividual(SEXP numberOfMarkersSEXP, SEXP numberOfAllelesSEXP, SEXP numberOfContributorsSEXP, SEXP numberOfKnownContributorsSEXP, SEXP knownProfilesSEXP, SEXP coverageSEXP, SEXP potentialParentsSEXP, SEXP markerImbalancesSEXP, SEXP convexMarkerImbalanceInterpolationSEXP, SEXP toleranceSEXP, SEXP thetaSEXP, SEXP alleleFrequenciesSEXP, SEXP levelsOfStutterRecursionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfMarkers(numberOfMarkersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type numberOfAlleles(numberOfAllelesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfContributors(numberOfContributorsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfKnownContributors(numberOfKnownContributorsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type knownProfiles(knownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector < Eigen::MatrixXd > >& >::type potentialParents(potentialParentsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type markerImbalances(markerImbalancesSEXP);
    Rcpp::traits::input_parameter< const double& >::type convexMarkerImbalanceInterpolation(convexMarkerImbalanceInterpolationSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alleleFrequencies(alleleFrequenciesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type levelsOfStutterRecursion(levelsOfStutterRecursionSEXP);
    rcpp_result_gen = Rcpp::wrap(setupIndividual(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, coverage, potentialParents, markerImbalances, convexMarkerImbalanceInterpolation, tolerance, theta, alleleFrequencies, levelsOfStutterRecursion));
    return rcpp_result_gen;
END_RCPP
}
// logPriorGenotypeProbability
Eigen::VectorXd logPriorGenotypeProbability(const Eigen::VectorXd& alleleFrequencies, const double& theta, const Eigen::MatrixXd& unknownProfiles, const Eigen::MatrixXd& knownProfiles, const std::size_t& numberOfMarkers, const Eigen::VectorXd& numberOfAlleles);
RcppExport SEXP _MPSMixtures_logPriorGenotypeProbability(SEXP alleleFrequenciesSEXP, SEXP thetaSEXP, SEXP unknownProfilesSEXP, SEXP knownProfilesSEXP, SEXP numberOfMarkersSEXP, SEXP numberOfAllelesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alleleFrequencies(alleleFrequenciesSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type unknownProfiles(unknownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type knownProfiles(knownProfilesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfMarkers(numberOfMarkersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type numberOfAlleles(numberOfAllelesSEXP);
    rcpp_result_gen = Rcpp::wrap(logPriorGenotypeProbability(alleleFrequencies, theta, unknownProfiles, knownProfiles, numberOfMarkers, numberOfAlleles));
    return rcpp_result_gen;
END_RCPP
}
// runningSinglePopulationEvolutionaryAlgorithm
Rcpp::List runningSinglePopulationEvolutionaryAlgorithm(const std::size_t& numberOfMarkers, const Eigen::VectorXd& numberOfAlleles, const std::size_t& numberOfContributors, const std::size_t& numberOfKnownContributors, const Eigen::MatrixXd& knownProfiles, const Eigen::MatrixXd& allKnownProfiles, const Eigen::VectorXd& coverage, const std::vector< std::vector < Eigen::MatrixXd > >& potentialParents, const Eigen::VectorXd& markerImbalances, const double& convexMarkerImbalanceInterpolation, const double& tolerance, const double& theta, const Eigen::VectorXd& alleleFrequencies, const std::size_t& populationSize, const std::size_t& numberOfIterations, const std::size_t& numberOfIterationsEqualMinMax, const std::size_t& numberOfFittestIndividuals, const int& parentSelectionWindowSize, const bool& allowParentSurvival, const double& crossoverProbability, const double& mutationProbabilityLowerLimit, const double& mutationDegreesOfFreedom, const Eigen::VectorXd mutationDecay, const std::size_t hillClimbingDirections, const std::size_t hillClimbingIterations, const std::size_t& seed, const bool& trace, const std::size_t& levelsOfStutterRecursion);
RcppExport SEXP _MPSMixtures_runningSinglePopulationEvolutionaryAlgorithm(SEXP numberOfMarkersSEXP, SEXP numberOfAllelesSEXP, SEXP numberOfContributorsSEXP, SEXP numberOfKnownContributorsSEXP, SEXP knownProfilesSEXP, SEXP allKnownProfilesSEXP, SEXP coverageSEXP, SEXP potentialParentsSEXP, SEXP markerImbalancesSEXP, SEXP convexMarkerImbalanceInterpolationSEXP, SEXP toleranceSEXP, SEXP thetaSEXP, SEXP alleleFrequenciesSEXP, SEXP populationSizeSEXP, SEXP numberOfIterationsSEXP, SEXP numberOfIterationsEqualMinMaxSEXP, SEXP numberOfFittestIndividualsSEXP, SEXP parentSelectionWindowSizeSEXP, SEXP allowParentSurvivalSEXP, SEXP crossoverProbabilitySEXP, SEXP mutationProbabilityLowerLimitSEXP, SEXP mutationDegreesOfFreedomSEXP, SEXP mutationDecaySEXP, SEXP hillClimbingDirectionsSEXP, SEXP hillClimbingIterationsSEXP, SEXP seedSEXP, SEXP traceSEXP, SEXP levelsOfStutterRecursionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfMarkers(numberOfMarkersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type numberOfAlleles(numberOfAllelesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfContributors(numberOfContributorsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfKnownContributors(numberOfKnownContributorsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type knownProfiles(knownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type allKnownProfiles(allKnownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector < Eigen::MatrixXd > >& >::type potentialParents(potentialParentsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type markerImbalances(markerImbalancesSEXP);
    Rcpp::traits::input_parameter< const double& >::type convexMarkerImbalanceInterpolation(convexMarkerImbalanceInterpolationSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alleleFrequencies(alleleFrequenciesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type populationSize(populationSizeSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfIterations(numberOfIterationsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfIterationsEqualMinMax(numberOfIterationsEqualMinMaxSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfFittestIndividuals(numberOfFittestIndividualsSEXP);
    Rcpp::traits::input_parameter< const int& >::type parentSelectionWindowSize(parentSelectionWindowSizeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type allowParentSurvival(allowParentSurvivalSEXP);
    Rcpp::traits::input_parameter< const double& >::type crossoverProbability(crossoverProbabilitySEXP);
    Rcpp::traits::input_parameter< const double& >::type mutationProbabilityLowerLimit(mutationProbabilityLowerLimitSEXP);
    Rcpp::traits::input_parameter< const double& >::type mutationDegreesOfFreedom(mutationDegreesOfFreedomSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type mutationDecay(mutationDecaySEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type hillClimbingDirections(hillClimbingDirectionsSEXP);
    Rcpp::traits::input_parameter< const std::size_t >::type hillClimbingIterations(hillClimbingIterationsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type levelsOfStutterRecursion(levelsOfStutterRecursionSEXP);
    rcpp_result_gen = Rcpp::wrap(runningSinglePopulationEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles, coverage, potentialParents, markerImbalances, convexMarkerImbalanceInterpolation, tolerance, theta, alleleFrequencies, populationSize, numberOfIterations, numberOfIterationsEqualMinMax, numberOfFittestIndividuals, parentSelectionWindowSize, allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, mutationDegreesOfFreedom, mutationDecay, hillClimbingDirections, hillClimbingIterations, seed, trace, levelsOfStutterRecursion));
    return rcpp_result_gen;
END_RCPP
}
// initialisingParallelEvolutionaryAlgorithm
Rcpp::List initialisingParallelEvolutionaryAlgorithm(const std::size_t& numberOfMarkers, const Eigen::VectorXd& numberOfAlleles, const std::size_t& numberOfContributors, const std::size_t& numberOfKnownContributors, const Eigen::MatrixXd& knownProfiles, const Eigen::MatrixXd& allKnownProfiles, const Eigen::VectorXd& coverage, const std::vector< std::vector < Eigen::MatrixXd > >& potentialParents, const Eigen::VectorXd& markerImbalances, const double& convexMarkerImbalanceInterpolation, const double& tolerance, const double& theta, const Eigen::VectorXd& alleleFrequencies, const std::size_t& populationSize, const std::size_t& seed, const std::size_t& levelsOfStutterRecursion);
RcppExport SEXP _MPSMixtures_initialisingParallelEvolutionaryAlgorithm(SEXP numberOfMarkersSEXP, SEXP numberOfAllelesSEXP, SEXP numberOfContributorsSEXP, SEXP numberOfKnownContributorsSEXP, SEXP knownProfilesSEXP, SEXP allKnownProfilesSEXP, SEXP coverageSEXP, SEXP potentialParentsSEXP, SEXP markerImbalancesSEXP, SEXP convexMarkerImbalanceInterpolationSEXP, SEXP toleranceSEXP, SEXP thetaSEXP, SEXP alleleFrequenciesSEXP, SEXP populationSizeSEXP, SEXP seedSEXP, SEXP levelsOfStutterRecursionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfMarkers(numberOfMarkersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type numberOfAlleles(numberOfAllelesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfContributors(numberOfContributorsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfKnownContributors(numberOfKnownContributorsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type knownProfiles(knownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type allKnownProfiles(allKnownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector < Eigen::MatrixXd > >& >::type potentialParents(potentialParentsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type markerImbalances(markerImbalancesSEXP);
    Rcpp::traits::input_parameter< const double& >::type convexMarkerImbalanceInterpolation(convexMarkerImbalanceInterpolationSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alleleFrequencies(alleleFrequenciesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type populationSize(populationSizeSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type levelsOfStutterRecursion(levelsOfStutterRecursionSEXP);
    rcpp_result_gen = Rcpp::wrap(initialisingParallelEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles, coverage, potentialParents, markerImbalances, convexMarkerImbalanceInterpolation, tolerance, theta, alleleFrequencies, populationSize, seed, levelsOfStutterRecursion));
    return rcpp_result_gen;
END_RCPP
}
// runningParallelEvolutionaryAlgorithm
Rcpp::List runningParallelEvolutionaryAlgorithm(const std::size_t& numberOfMarkers, const Eigen::VectorXd& numberOfAlleles, const std::size_t& numberOfContributors, const std::size_t& numberOfKnownContributors, const Eigen::MatrixXd& knownProfiles, const Eigen::MatrixXd& allKnownProfiles, const Eigen::VectorXd& coverage, const std::vector< std::vector < Eigen::MatrixXd > >& potentialParents, const Eigen::VectorXd& markerImbalances, const double& convexMarkerImbalanceInterpolation, const double& tolerance, const double& theta, const Eigen::VectorXd& alleleFrequencies, const std::size_t& numberOfIterations, const std::size_t& numberOfIterationsEqualMinMax, const std::size_t& numberOfFittestIndividuals, const int& parentSelectionWindowSize, const bool allowParentSurvival, const double& crossoverProbability, const double& mutationProbabilityLowerLimit, const double& mutationDegreesOfFreedom, const Eigen::VectorXd& mutationDecay, const std::size_t& hillClimbingDirections, const std::size_t& hillClimbingIterations, const std::size_t& seed, const bool& trace, const Eigen::MatrixXd encodedPopulationList, const Eigen::MatrixXd sampleParametersList, const Eigen::MatrixXd noiseParametersList, const Eigen::MatrixXd mixtureParametersList, const Eigen::MatrixXd markerParametersList, const Eigen::VectorXd fitnessList, const std::size_t& levelsOfStutterRecursion);
RcppExport SEXP _MPSMixtures_runningParallelEvolutionaryAlgorithm(SEXP numberOfMarkersSEXP, SEXP numberOfAllelesSEXP, SEXP numberOfContributorsSEXP, SEXP numberOfKnownContributorsSEXP, SEXP knownProfilesSEXP, SEXP allKnownProfilesSEXP, SEXP coverageSEXP, SEXP potentialParentsSEXP, SEXP markerImbalancesSEXP, SEXP convexMarkerImbalanceInterpolationSEXP, SEXP toleranceSEXP, SEXP thetaSEXP, SEXP alleleFrequenciesSEXP, SEXP numberOfIterationsSEXP, SEXP numberOfIterationsEqualMinMaxSEXP, SEXP numberOfFittestIndividualsSEXP, SEXP parentSelectionWindowSizeSEXP, SEXP allowParentSurvivalSEXP, SEXP crossoverProbabilitySEXP, SEXP mutationProbabilityLowerLimitSEXP, SEXP mutationDegreesOfFreedomSEXP, SEXP mutationDecaySEXP, SEXP hillClimbingDirectionsSEXP, SEXP hillClimbingIterationsSEXP, SEXP seedSEXP, SEXP traceSEXP, SEXP encodedPopulationListSEXP, SEXP sampleParametersListSEXP, SEXP noiseParametersListSEXP, SEXP mixtureParametersListSEXP, SEXP markerParametersListSEXP, SEXP fitnessListSEXP, SEXP levelsOfStutterRecursionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfMarkers(numberOfMarkersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type numberOfAlleles(numberOfAllelesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfContributors(numberOfContributorsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfKnownContributors(numberOfKnownContributorsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type knownProfiles(knownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type allKnownProfiles(allKnownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector < Eigen::MatrixXd > >& >::type potentialParents(potentialParentsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type markerImbalances(markerImbalancesSEXP);
    Rcpp::traits::input_parameter< const double& >::type convexMarkerImbalanceInterpolation(convexMarkerImbalanceInterpolationSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alleleFrequencies(alleleFrequenciesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfIterations(numberOfIterationsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfIterationsEqualMinMax(numberOfIterationsEqualMinMaxSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfFittestIndividuals(numberOfFittestIndividualsSEXP);
    Rcpp::traits::input_parameter< const int& >::type parentSelectionWindowSize(parentSelectionWindowSizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type allowParentSurvival(allowParentSurvivalSEXP);
    Rcpp::traits::input_parameter< const double& >::type crossoverProbability(crossoverProbabilitySEXP);
    Rcpp::traits::input_parameter< const double& >::type mutationProbabilityLowerLimit(mutationProbabilityLowerLimitSEXP);
    Rcpp::traits::input_parameter< const double& >::type mutationDegreesOfFreedom(mutationDegreesOfFreedomSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mutationDecay(mutationDecaySEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type hillClimbingDirections(hillClimbingDirectionsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type hillClimbingIterations(hillClimbingIterationsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type encodedPopulationList(encodedPopulationListSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type sampleParametersList(sampleParametersListSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type noiseParametersList(noiseParametersListSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type mixtureParametersList(mixtureParametersListSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type markerParametersList(markerParametersListSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type fitnessList(fitnessListSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type levelsOfStutterRecursion(levelsOfStutterRecursionSEXP);
    rcpp_result_gen = Rcpp::wrap(runningParallelEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles, coverage, potentialParents, markerImbalances, convexMarkerImbalanceInterpolation, tolerance, theta, alleleFrequencies, numberOfIterations, numberOfIterationsEqualMinMax, numberOfFittestIndividuals, parentSelectionWindowSize, allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, mutationDegreesOfFreedom, mutationDecay, hillClimbingDirections, hillClimbingIterations, seed, trace, encodedPopulationList, sampleParametersList, noiseParametersList, mixtureParametersList, markerParametersList, fitnessList, levelsOfStutterRecursion));
    return rcpp_result_gen;
END_RCPP
}
// oneStepApproximationCpp
Rcpp::List oneStepApproximationCpp(const Eigen::VectorXd& encodedProfiles, const Eigen::VectorXd& sampleParameters, const Eigen::VectorXd& noiseParameters, const Eigen::VectorXd& mixtureParameters, const Eigen::VectorXd& coverage, const Eigen::VectorXd markerImbalances, const std::vector< std::vector< Eigen::MatrixXd > > potentialParents, const Eigen::MatrixXd& knownProfiles, const Eigen::MatrixXd& allKnownProfiles, const Eigen::VectorXd& alleleFrequencies, const double& theta, const std::size_t& numberOfContributors, const std::size_t& numberOfMarkers, const Eigen::VectorXd& numberOfAlleles, const std::size_t& levelsOfStutterRecursion);
RcppExport SEXP _MPSMixtures_oneStepApproximationCpp(SEXP encodedProfilesSEXP, SEXP sampleParametersSEXP, SEXP noiseParametersSEXP, SEXP mixtureParametersSEXP, SEXP coverageSEXP, SEXP markerImbalancesSEXP, SEXP potentialParentsSEXP, SEXP knownProfilesSEXP, SEXP allKnownProfilesSEXP, SEXP alleleFrequenciesSEXP, SEXP thetaSEXP, SEXP numberOfContributorsSEXP, SEXP numberOfMarkersSEXP, SEXP numberOfAllelesSEXP, SEXP levelsOfStutterRecursionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type encodedProfiles(encodedProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type sampleParameters(sampleParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type noiseParameters(noiseParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mixtureParameters(mixtureParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type markerImbalances(markerImbalancesSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< Eigen::MatrixXd > > >::type potentialParents(potentialParentsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type knownProfiles(knownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type allKnownProfiles(allKnownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alleleFrequencies(alleleFrequenciesSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfContributors(numberOfContributorsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfMarkers(numberOfMarkersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type numberOfAlleles(numberOfAllelesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type levelsOfStutterRecursion(levelsOfStutterRecursionSEXP);
    rcpp_result_gen = Rcpp::wrap(oneStepApproximationCpp(encodedProfiles, sampleParameters, noiseParameters, mixtureParameters, coverage, markerImbalances, potentialParents, knownProfiles, allKnownProfiles, alleleFrequencies, theta, numberOfContributors, numberOfMarkers, numberOfAlleles, levelsOfStutterRecursion));
    return rcpp_result_gen;
END_RCPP
}
// EAApproximationCpp
Rcpp::List EAApproximationCpp(const std::vector<Eigen::VectorXd>& encodedProfiles, const Eigen::VectorXd& sampleParameters, const Eigen::VectorXd& noiseParameters, const Eigen::VectorXd& mixtureParameters, const Eigen::VectorXd& coverage, const Eigen::VectorXd markerImbalances, const std::vector< std::vector< Eigen::MatrixXd > > potentialParents, const Eigen::MatrixXd& knownProfiles, const Eigen::MatrixXd& allKnownProfiles, const Eigen::VectorXd& alleleFrequencies, const double& theta, const std::size_t& numberOfContributors, const std::size_t& numberOfMarkers, const Eigen::VectorXd& numberOfAlleles, const std::size_t& levelsOfStutterRecursion);
RcppExport SEXP _MPSMixtures_EAApproximationCpp(SEXP encodedProfilesSEXP, SEXP sampleParametersSEXP, SEXP noiseParametersSEXP, SEXP mixtureParametersSEXP, SEXP coverageSEXP, SEXP markerImbalancesSEXP, SEXP potentialParentsSEXP, SEXP knownProfilesSEXP, SEXP allKnownProfilesSEXP, SEXP alleleFrequenciesSEXP, SEXP thetaSEXP, SEXP numberOfContributorsSEXP, SEXP numberOfMarkersSEXP, SEXP numberOfAllelesSEXP, SEXP levelsOfStutterRecursionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<Eigen::VectorXd>& >::type encodedProfiles(encodedProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type sampleParameters(sampleParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type noiseParameters(noiseParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mixtureParameters(mixtureParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type markerImbalances(markerImbalancesSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< Eigen::MatrixXd > > >::type potentialParents(potentialParentsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type knownProfiles(knownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type allKnownProfiles(allKnownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alleleFrequencies(alleleFrequenciesSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfContributors(numberOfContributorsSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfMarkers(numberOfMarkersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type numberOfAlleles(numberOfAllelesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type levelsOfStutterRecursion(levelsOfStutterRecursionSEXP);
    rcpp_result_gen = Rcpp::wrap(EAApproximationCpp(encodedProfiles, sampleParameters, noiseParameters, mixtureParameters, coverage, markerImbalances, potentialParents, knownProfiles, allKnownProfiles, alleleFrequencies, theta, numberOfContributors, numberOfMarkers, numberOfAlleles, levelsOfStutterRecursion));
    return rcpp_result_gen;
END_RCPP
}
// samplePosteriorGenotypesGuidedCpp
Rcpp::List samplePosteriorGenotypesGuidedCpp(const Eigen::VectorXd& encodedProfiles, const Eigen::VectorXd& sampleParameters, const Eigen::VectorXd& noiseParameters, const Eigen::VectorXd& mixtureParameters, const Eigen::VectorXd markerParameters, const Eigen::VectorXd& coverage, const std::vector< std::vector< Eigen::MatrixXd > > potentialParents, const Eigen::MatrixXd& knownProfiles, const Eigen::MatrixXd& allKnownProfiles, const Eigen::VectorXd& alleleFrequencies, const double& theta, const std::size_t& numberOfContributors, const Eigen::VectorXd& numberOfAlleles, const std::size_t& levelsOfStutterRecursion, const std::size_t& numberOfSimulations, const bool& suggestionBool, const std::size_t& seed);
RcppExport SEXP _MPSMixtures_samplePosteriorGenotypesGuidedCpp(SEXP encodedProfilesSEXP, SEXP sampleParametersSEXP, SEXP noiseParametersSEXP, SEXP mixtureParametersSEXP, SEXP markerParametersSEXP, SEXP coverageSEXP, SEXP potentialParentsSEXP, SEXP knownProfilesSEXP, SEXP allKnownProfilesSEXP, SEXP alleleFrequenciesSEXP, SEXP thetaSEXP, SEXP numberOfContributorsSEXP, SEXP numberOfAllelesSEXP, SEXP levelsOfStutterRecursionSEXP, SEXP numberOfSimulationsSEXP, SEXP suggestionBoolSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type encodedProfiles(encodedProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type sampleParameters(sampleParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type noiseParameters(noiseParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mixtureParameters(mixtureParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type markerParameters(markerParametersSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type coverage(coverageSEXP);
    Rcpp::traits::input_parameter< const std::vector< std::vector< Eigen::MatrixXd > > >::type potentialParents(potentialParentsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type knownProfiles(knownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type allKnownProfiles(allKnownProfilesSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type alleleFrequencies(alleleFrequenciesSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfContributors(numberOfContributorsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type numberOfAlleles(numberOfAllelesSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type levelsOfStutterRecursion(levelsOfStutterRecursionSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type numberOfSimulations(numberOfSimulationsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type suggestionBool(suggestionBoolSEXP);
    Rcpp::traits::input_parameter< const std::size_t& >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(samplePosteriorGenotypesGuidedCpp(encodedProfiles, sampleParameters, noiseParameters, mixtureParameters, markerParameters, coverage, potentialParents, knownProfiles, allKnownProfiles, alleleFrequencies, theta, numberOfContributors, numberOfAlleles, levelsOfStutterRecursion, numberOfSimulations, suggestionBool, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MPSMixtures_partialSumEigen", (DL_FUNC) &_MPSMixtures_partialSumEigen, 1},
    {"_MPSMixtures_setupIndividual", (DL_FUNC) &_MPSMixtures_setupIndividual, 13},
    {"_MPSMixtures_logPriorGenotypeProbability", (DL_FUNC) &_MPSMixtures_logPriorGenotypeProbability, 6},
    {"_MPSMixtures_runningSinglePopulationEvolutionaryAlgorithm", (DL_FUNC) &_MPSMixtures_runningSinglePopulationEvolutionaryAlgorithm, 28},
    {"_MPSMixtures_initialisingParallelEvolutionaryAlgorithm", (DL_FUNC) &_MPSMixtures_initialisingParallelEvolutionaryAlgorithm, 16},
    {"_MPSMixtures_runningParallelEvolutionaryAlgorithm", (DL_FUNC) &_MPSMixtures_runningParallelEvolutionaryAlgorithm, 33},
    {"_MPSMixtures_oneStepApproximationCpp", (DL_FUNC) &_MPSMixtures_oneStepApproximationCpp, 15},
    {"_MPSMixtures_EAApproximationCpp", (DL_FUNC) &_MPSMixtures_EAApproximationCpp, 15},
    {"_MPSMixtures_samplePosteriorGenotypesGuidedCpp", (DL_FUNC) &_MPSMixtures_samplePosteriorGenotypesGuidedCpp, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_MPSMixtures(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
