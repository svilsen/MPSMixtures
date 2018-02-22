// estimatePoissonGammaParameters.hpp
#ifndef estimatePoissonGammaParameters
#define estimatePoissonGammaParameters

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include "experimentalSetup.hpp"

// Allele coverage
class EstimatePoissonGammaAlleleParameters
{
    public:
        std::size_t NumberOfContributors;

        Eigen::VectorXd Coverage;
        Eigen::VectorXd MarkerImbalances;
        Eigen::MatrixXd ExpectedContributionMatrix;
        Eigen::VectorXd NumberOfAlleles;
        Eigen::VectorXd AlleleCount;

        Eigen::VectorXd SampleParameters;
        Eigen::VectorXd MixtureParameters;
        Eigen::VectorXd MarkerImbalancesParameters;

        double ConvexMarkerImbalanceInterpolation;

        double LogLikelihood;
        double ToleranceFRelative;
        double ToleranceEqualityConstraints;
        unsigned int MaximumNumberOfIterations;

        std::size_t Counter;

        EstimatePoissonGammaAlleleParameters(Eigen::VectorXd coverage, Eigen::MatrixXd expectedContributionMatrix,
                                             Eigen::VectorXd markerImbalances, Eigen::VectorXd numberOfAlleles, Eigen::VectorXd alleleCount,
                                             double convexMarkerImbalanceInterpolation, double tolerance);
        void initialiseParameters();
};

double logLikelihoodAlleleCoverage(const std::vector<double> &x, std::vector<double> &grad, void *data);
double linearConstraintsLoglikelihoodAlleleCoverage(const std::vector<double>&x, std::vector<double> &grad, void *data);

void estimateParametersAlleleCoverage(EstimatePoissonGammaAlleleParameters &EPGA);

// Noise coverage
class EstimatePoissonGammaNoiseParameters
{
    public:
        Eigen::VectorXd Coverage;
        Eigen::VectorXd NoiseParameters;

        double LogLikelihood;
        double ToleranceFRelative;
        unsigned int MaximumNumberOfIterations;

        std::size_t Counter;

        EstimatePoissonGammaNoiseParameters(Eigen::VectorXd coverage, double tolerance);
        void initialiseParameters();
};

double logLikelihoodNoiseCoverage(const std::vector<double> &x, std::vector<double> &grad, void *data);

void estimateParametersNoiseCoverage(EstimatePoissonGammaNoiseParameters &EPGN);


#endif
