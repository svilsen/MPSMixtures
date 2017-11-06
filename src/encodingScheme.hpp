// encodingScheme
#ifndef encodingScheme
#define encodingScheme

#include <RcppEigen.h>

Eigen::MatrixXd decoding(Eigen::VectorXd individual, Eigen::VectorXd numberOfAlleles, std::size_t numberOfMarkers, std::size_t numberOfUnknownContributors);

#endif
