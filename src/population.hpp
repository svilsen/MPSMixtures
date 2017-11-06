// population.hpp
#ifndef population
#define population

#include <RcppEigen.h>
#include <vector>

#include "individual.hpp"

class Population
{
    public:
        std::vector<Individual> Individuals;
        // Eigen::MatrixXd EncodedProfiles;
        Eigen::VectorXd Fitness;

        Population();
        Population(std::vector<Individual> I);

        Rcpp::List ReturnRcppList();
        Rcpp::List ReturnCompressedRcppList();
};

#endif
