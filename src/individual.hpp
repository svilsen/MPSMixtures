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
        std::vector<Eigen::MatrixXd> ReducedExpectedContributionMatrix;
        std::vector<Eigen::VectorXd> ReducedAlleleIndex;
        std::vector<Eigen::VectorXd> ReducedNoiseIndex;

        Eigen::VectorXd SampleParameters;
        Eigen::VectorXd MixtureParameters;
        Eigen::VectorXd MarkerImbalanceParameters;

        Eigen::VectorXd NoiseParameters;

        // Likelihoods and fitness
        Eigen::VectorXd LogLikelihoodAlleleMarker;
        double LogLikelihoodAllele;

        Eigen::VectorXd LogLikelihoodNoiseMarker;
        double LogLikelihoodNoise;

        Eigen::VectorXd LogPriorGenotypeProbabilityMarker;
        double LogPriorGenotypeProbability;

        double Fitness;

        // Constructors
        Individual();

        Individual(const ExperimentalSetup & ES);

        Individual(const Eigen::VectorXd & encodedGenotype, const ExperimentalSetup & ES);

        Individual(const Eigen::VectorXd & encodedGenotype, const Eigen::VectorXd & sampleParameters, const Eigen::VectorXd & noiseParameters,
                   const Eigen::VectorXd & mixtureParameters, const Eigen::VectorXd & markerParameters, const double & fitness);

        Individual(const Eigen::VectorXd & encodedGenotype, const Eigen::VectorXd & sampleParameters, const Eigen::VectorXd & noiseParameters,
                   const Eigen::VectorXd & mixtureParameters, const Eigen::VectorXd & markerParameters, const ExperimentalSetup & ES);

        Individual(const Eigen::VectorXd & encodedGenotype, const std::vector<Eigen::MatrixXd> reducedExpectedContributionMatrix,
                   const ExperimentalSetup & ES);

        // ECM, Genotype, noise, etc. functions
        void CreateReducedElements(const ExperimentalSetup & ES);

        Eigen::MatrixXd GenerateExpectedContributionProfile(const ExperimentalSetup & ES, const Eigen::MatrixXd & decodedProfile);
        Eigen::VectorXd GenerateNoiseProfile(const ExperimentalSetup & ES, const Eigen::MatrixXd & expectedContributionProfile);

        //
        void EstimateParameters(const ExperimentalSetup & ES, Eigen::MatrixXd decodedProfile, Eigen::MatrixXd expectedContributionProfile,
                                Eigen::VectorXd noiseProfile);

        Eigen::VectorXd CalculateResiduals(const ExperimentalSetup & ES, const Eigen::MatrixXd & expectedContributionProfile, const Eigen::VectorXd & noiseProfile);
        void CalculateFitness(const ExperimentalSetup & ES, Eigen::MatrixXd decodedProfile, Eigen::MatrixXd expectedContributionProfile,
                              Eigen::VectorXd noiseProfile);

        // R-wrappers
        Rcpp::List ReturnRcppList(const ExperimentalSetup & ES);
        Rcpp::List ReturnRcppListSimplified();
};

double ParentStutterContribution(const std::size_t & currentAllele, const std::size_t & stutterRecursion, const double & levelOfRecursion,
                                 const Eigen::VectorXd & decodedProfile_mu, const std::vector<Eigen::MatrixXd> & potentialParents_m, const double & numberOfAlleles_m);

#endif

