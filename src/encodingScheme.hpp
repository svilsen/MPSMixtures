// encodingScheme
#ifndef encodingScheme
#define encodingScheme

#include <RcppEigen.h>

Eigen::MatrixXd decoding(const Eigen::VectorXd & individual, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfMarkers,
                         const std::size_t & numberOfUnknownContributors);

#endif
