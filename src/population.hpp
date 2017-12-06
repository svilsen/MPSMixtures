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
        // Eigen::MatrixXd EncodedProfiles;
        Eigen::VectorXd Fitness;

        Population();
        Population(std::vector<Individual> I);
        Population(Eigen::MatrixXd encodedProfilesList, Eigen::MatrixXd sampleParametersList, Eigen::MatrixXd noiseParametersList,
                   Eigen::MatrixXd mixtureParametersList, ExperimentalSetup ES);

        Rcpp::List ReturnRcppList();
        Rcpp::List ReturnCompressedRcppList();
};

#endif
