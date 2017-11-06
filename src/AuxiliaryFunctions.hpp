// AuxiliaryFunctions.hpp
#ifndef AuxiliaryFunctions
#define AuxiliaryFunctions

#include <RcppEigen.h>

Eigen::VectorXd partialSumEigen(const Eigen::VectorXd & x);

Eigen::MatrixXd bindColumns(const Eigen::MatrixXd & A, const Eigen::MatrixXd & B);

Eigen::VectorXd STDEigen(std::vector<double> &x);
std::vector<double> EigenSTD(Eigen::VectorXd &x);

Eigen::Vector2i nonZeroElementsOfMarker(const Eigen::VectorXd & decodedProfile_m);

std::vector<int> sortedIndex(const Eigen::VectorXd & x);

#endif
