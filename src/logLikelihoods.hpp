// logLikelihoods.hpp
#ifndef logLikelihoods
#define logLikelihoods

#include <RcppEigen.h>

double logPoissonGammaDistribution(const double & x, const double & mean, const double & overdispersion);
double devianceResidualPoissonGammaDistribution(const double & x, const double & mean, const double & overdispersion);

double devianceResidualPG1(const double & x, const double & mean, const double & dispersion);

double logMultinomialCoefficient(const int & totalCounts, const Eigen::VectorXd & counts);
double logDirichletMultinomial(const Eigen::VectorXd & counts, const Eigen::VectorXd & alleleFrequencies);
double logDirichletMultinomialTheta(const Eigen::VectorXd & counts, const double & theta, const Eigen::VectorXd & alleleFrequencies);

Eigen::VectorXd logLikelihoodAlleleCoverage(const Eigen::VectorXd & coverage, const Eigen::MatrixXd & expectedContributionMatrix,
                                            const Eigen::VectorXd & sampleParameters, const Eigen::VectorXd & mixtureProportions,
                                            const Eigen::VectorXd & numberOfAlleles, const Eigen::VectorXd & markerImbalances);

Eigen::VectorXd logLikelihoodNoiseCoverage(const Eigen::VectorXd & coverage, const Eigen::VectorXd & noiseProfile, const double & noiseLevel, const double & noiseDispersion, const Eigen::VectorXd & numberOfAlleles);

Eigen::VectorXd logPriorGenotypeProbability(const Eigen::VectorXd & alleleFrequencies, const double & theta, const Eigen::MatrixXd & unknownProfiles, const Eigen::MatrixXd & knownProfiles, const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles);

#endif
