// estimatePoissonGammaAlleleCoverageParameters.hpp
#ifndef estimatePoissonGammaAlleleCoverageParameters
#define estimatePoissonGammaAlleleCoverageParameters

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

std::vector<Eigen::VectorXd> estimateParametersAlleleCoverage(const Eigen::VectorXd & coverage, const Eigen::MatrixXd & expectedContributionMatrix,
                                                              const Eigen::VectorXd & markerImbalances, const Eigen::VectorXd & numberOfAlleles,
                                                              const Eigen::VectorXd & alleleCount, const double & convexMarkerImbalanceInterpolation,
                                                              const double & tolerance);

Eigen::VectorXd estimateParametersNoiseCoverage(const Eigen::VectorXd & coverage);

#endif
