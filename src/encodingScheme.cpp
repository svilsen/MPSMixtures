#include <Rcpp.h>
#include <RcppEigen.h>

#include "AuxiliaryFunctions.hpp"

Eigen::MatrixXd decoding(Eigen::VectorXd individual, Eigen::VectorXd numberOfAlleles, std::size_t numberOfMarkers, std::size_t numberOfUnknownContributors)
{
    const std::size_t N = individual.size();
    const Eigen::VectorXd partialSumAlleles = partialSumEigen(numberOfAlleles);

    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(numberOfAlleles.sum(), numberOfUnknownContributors);
    for (std::size_t i = 0; i < numberOfMarkers; i++)
    {
        for (std::size_t j = 0; j < 2*numberOfUnknownContributors; j++)
        {
            const std::size_t k = 2*numberOfUnknownContributors*i + j;

            std::size_t column = std::floor(j/2);
            std::size_t row = individual[k] + partialSumAlleles[i];
            res(row, column) += 1;
        }
    }

    return res;
}
