// population.hpp
#ifndef population
#define population

#include <RcppEigen.h>
#include <vector>

#include "experimentalSetup.hpp"
#include "individual.hpp"

class Population
{
    public:
        std::vector<Individual> Individuals;
        Eigen::VectorXd Fitness;

        Population();
        Population(std::vector<Individual> I);
        Population(const Eigen::MatrixXd & encodedProfilesList, const Eigen::MatrixXd & sampleParametersList, const Eigen::MatrixXd & noiseParametersList,
                   const Eigen::MatrixXd & mixtureParametersList, const Eigen::MatrixXd & markerParametersList, const Eigen::VectorXd & fitnessList);

        Rcpp::List ReturnRcppList(const ExperimentalSetup & ES);
        Rcpp::List ReturnCompressedRcppList();
};

#endif
