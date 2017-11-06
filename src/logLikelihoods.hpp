// logLikelihoods.hpp
#ifndef logLikelihoods
#define logLikelihoods

#include <RcppEigen.h>

double logPoissonGammaDistribution(double & x, double & mean, double & overdispersion);

double logMultinomialCoefficient(int totalCounts, Eigen::VectorXd counts);
double logDirichletMultinomial(Eigen::VectorXd counts, Eigen::VectorXd alleleFrequencies);
double logDirichletMultinomialTheta(Eigen::VectorXd counts, double theta, Eigen::VectorXd alleleFrequencies);

double logLikelihoodAlleleCoverage(Eigen::VectorXd coverage, Eigen::MatrixXd expectedContributionMatrix,
                                   Eigen::VectorXd sampleParameters, Eigen::VectorXd mixtureProportions,
                                   Eigen::VectorXd markerImbalances);

double logLikelihoodNoiseCoverage(Eigen::VectorXd & coverage, Eigen::VectorXd & noiseProfile, double noiseLevel, double noiseDispersion);

double logPriorGenotypeProbability(Eigen::VectorXd & alleleFrequencies, double theta, Eigen::MatrixXd & unknownProfiles, Eigen::MatrixXd & knownProfiles, std::size_t numberOfMarkers, Eigen::VectorXd numberOfAlleles);

#endif
