#include <Rcpp.h>
#include <RcppEigen.h>

#include "population.hpp"
#include "individual.hpp"

Population::Population() { }

Population::Population(std::vector<Individual> I)
{
    Individuals = I;

    // std::size_t N = Individuals[0].EncodedProfile.size();
    std::size_t M = Individuals.size();

    // Eigen::MatrixXd encodedProfileMatrix = Eigen::MatrixXd::Zero(N, M);
    Eigen::VectorXd fitnessVector = Eigen::VectorXd::Zero(M);
    for (std::size_t m = 0; m < M; m++)
    {
        // encodedProfileMatrix.col(m) = Individuals[m].EncodedProfile;
        fitnessVector[m] = Individuals[m].Fitness;
    }

    // EncodedProfiles = encodedProfileMatrix;
    Fitness = fitnessVector;
}

Population::Population(Eigen::MatrixXd encodedProfilesList, Eigen::MatrixXd sampleParametersList, Eigen::MatrixXd noiseParametersList,
                       Eigen::MatrixXd mixtureParametersList, ExperimentalSetup ES)
{
    const std::size_t populationSize = encodedProfilesList.cols();

    std::vector< Individual > currentIndividuals(populationSize);
    for (std::size_t i = 0; i < populationSize; i++)
    {
        Eigen::VectorXd encodedProfile_i = encodedProfilesList.col(i);
        Eigen::VectorXd sampleParameters_i = sampleParametersList.col(i);
        Eigen::VectorXd noiseParameters_i = noiseParametersList.col(i);
        Eigen::VectorXd mixtureParameters_i = mixtureParametersList.col(i);
        Individual I(encodedProfile_i, sampleParameters_i, noiseParameters_i, mixtureParameters_i, ES);
        currentIndividuals[i] = I;
    }

    Individuals = currentIndividuals;

    Eigen::VectorXd fitnessVector = Eigen::VectorXd::Zero(populationSize);
    for (std::size_t m = 0; m < populationSize; m++)
    {
        fitnessVector[m] = Individuals[m].Fitness;
    }

    Fitness = fitnessVector;
}


Rcpp::List Population::ReturnRcppList()
{
    std::size_t populationSize = Individuals.size();

    Rcpp::List RL(populationSize);
    for (std::size_t n = 0; n < populationSize; n++)
    {
        RL[n] = Individuals[n].ReturnRcppList();
    }

    return RL;
}

Rcpp::List Population::ReturnCompressedRcppList()
{
    std::size_t populationSize = Individuals.size();

    Eigen::MatrixXd encodedIndividuals = Eigen::MatrixXd::Zero(Individuals[0].EncodedProfile.size(), populationSize);
    Eigen::MatrixXd sampleParameters = Eigen::MatrixXd::Zero(Individuals[0].SampleParameters.size(), populationSize);
    Eigen::MatrixXd noiseParameters = Eigen::MatrixXd::Zero(Individuals[0].NoiseParameters.size(), populationSize);
    Eigen::MatrixXd mixtureParameters = Eigen::MatrixXd::Zero(Individuals[0].MixtureParameters.size(), populationSize);
    for (std::size_t p = 0; p < populationSize; p++)
    {
        Individual I_m = Individuals[p];
        encodedIndividuals.col(p) = I_m.EncodedProfile;
        sampleParameters.col(p) = I_m.SampleParameters;
        noiseParameters.col(p) = I_m.NoiseParameters;
        mixtureParameters.col(p) = I_m.MixtureParameters;
    }

    return Rcpp::List::create(Rcpp::Named("EncodedProfiles") = encodedIndividuals,
                              Rcpp::Named("SampleParameters") = sampleParameters,
                              Rcpp::Named("NoiseParameters") = noiseParameters,
                              Rcpp::Named("MixtureParameters") = mixtureParameters,
                              Rcpp::Named("Fitness") = Fitness);
}
