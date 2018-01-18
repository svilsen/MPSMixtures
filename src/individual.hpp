// individual.hpp
#ifndef individual
#define individual

#include <RcppEigen.h>

#include "experimentalSetup.hpp"
#include "estimatePoissonGammaParameters.hpp"

class Individual
{
    public:
        // Objects
        Eigen::VectorXd EncodedProfile;

        Eigen::VectorXd SampleParameters;
        Eigen::VectorXd MixtureParameters;
        Eigen::VectorXd NoiseParameters;


        Eigen::VectorXd LogLikelihoodAlleleMarker;
        double LogLikelihoodAllele;
        double LogLikelihoodNoise;
        double LogPriorGenotypeProbability;

        double Fitness;

        Individual();
        Individual(const ExperimentalSetup & ES);
        Individual(const Eigen::VectorXd & encodedGenotype, const ExperimentalSetup & ES);
        Individual(const Eigen::VectorXd & encodedGenotype, const Eigen::VectorXd & sampleParameters, const Eigen::VectorXd & noiseParameters,
                   const Eigen::VectorXd & mixtureParameters, const double & fitness);

        Eigen::MatrixXd GenerateExpectedContributionProfile(const ExperimentalSetup & ES, const Eigen::MatrixXd & decodedProfile);
        Eigen::VectorXd GenerateNoiseProfile(const ExperimentalSetup & ES, const Eigen::MatrixXd & expectedContributionProfile);

        void EstimateParameters(const ExperimentalSetup & ES);
        Eigen::VectorXd CalculateResiduals(const ExperimentalSetup & ES, const Eigen::MatrixXd & expectedContributionProfile, const Eigen::VectorXd & noiseProfile);
        void CalculateFitness(const ExperimentalSetup & ES);

        Rcpp::List ReturnRcppList(const ExperimentalSetup & ES);
        Rcpp::List ReturnRcppListSimplified();
};

#endif

