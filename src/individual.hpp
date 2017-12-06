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
        Eigen::MatrixXd DecodedProfile;
        Eigen::MatrixXd ExpectedContributionProfile;
        Eigen::VectorXd NoiseProfile;

        Eigen::VectorXd SampleParameters;
        Eigen::VectorXd MixtureParameters;
        Eigen::VectorXd NoiseParameters;

        Eigen::VectorXd DevianceResiduals;

        Eigen::VectorXd LogLikelihoodAlleleMarker;
        double LogLikelihoodAllele;
        double LogLikelihoodNoise;
        double LogPriorGenotypeProbability;

        double Fitness;

        // Functions
        Individual();
        Individual(ExperimentalSetup ES);
        Individual(Eigen::VectorXd encodedGenotype, ExperimentalSetup ES);
        Individual(Eigen::VectorXd encodedGenotype, Eigen::VectorXd sampleParameters, Eigen::VectorXd noiseParameters,
                   Eigen::VectorXd mixtureParameters, ExperimentalSetup ES);

        Eigen::MatrixXd GenerateExpectedContributionProfile(ExperimentalSetup ES);
        Eigen::VectorXd GenerateNoiseProfile(ExperimentalSetup ES);

        void EstimateParameters(ExperimentalSetup ES);
        void CalculateResiduals(ExperimentalSetup ES);
        void CalculateFitness(ExperimentalSetup ES);

        Rcpp::List ReturnRcppList();
        Rcpp::List ReturnRcppListSimplified();
};

#endif

