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
        std::size_t NumberOfMarkers;

        Eigen::VectorXd Coverage;
        Eigen::VectorXd MarkerImbalances;
        Eigen::VectorXd PartialSumAlleles;

        std::vector<Eigen::MatrixXd> ExpectedContributionMatrix;
        std::vector<Eigen::VectorXd> AlleleIndex;

        Eigen::VectorXd SampleParameters;
        Eigen::VectorXd MixtureParameters;
        Eigen::VectorXd MarkerImbalancesParameters;

        double ConvexMarkerImbalanceInterpolation;

        double LogLikelihood;
        Eigen::VectorXd Tolerance;
        unsigned int MaximumNumberOfIterations;

        std::size_t Counter;

        EstimatePoissonGammaAlleleParameters(const Eigen::VectorXd & coverage, const std::vector<Eigen::MatrixXd> & expectedContributionMatrix,
                                             const std::vector<Eigen::VectorXd> & alleleIndex,
                                             const Eigen::VectorXd & markerImbalances, const Eigen::VectorXd & partialSumAlleles,
                                             const double & convexMarkerImbalanceInterpolation, const Eigen::VectorXd & tolerance);

        void initialiseParameters();
};

double logLikelihoodAlleleCoverageNLopt(const std::vector<double> &x, std::vector<double> &grad, void *data);
double linearConstraintsLoglikelihoodAlleleCoverage(const std::vector<double>&x, std::vector<double> &grad, void *data);

void estimateParametersAlleleCoverage(EstimatePoissonGammaAlleleParameters &EPGA);

// Noise coverage
class EstimatePoissonGammaNoiseParameters
{
    public:
        std::size_t NumberOfMarkers;

        Eigen::VectorXd Coverage;
        std::vector<Eigen::VectorXd> NoiseIndex;
        Eigen::VectorXd PartialSumAlleles;

        Eigen::VectorXd NoiseParameters;
        double VarianceUpperLimit;

        double LogLikelihood;
        Eigen::VectorXd Tolerance;
        unsigned int MaximumNumberOfIterations;

        std::size_t Counter;

        EstimatePoissonGammaNoiseParameters(const Eigen::VectorXd & coverage, const std::vector<Eigen::VectorXd> & noiseIndex,
                                            const Eigen::VectorXd & partialSumAlleles, const Eigen::VectorXd & tolerance,
                                            const double & varianceUpperLimit);
        void initialiseParameters();
};

double logLikelihoodNoiseCoverageNLopt(const std::vector<double> &x, std::vector<double> &grad, void *data);
void estimateParametersNoiseCoverage(EstimatePoissonGammaNoiseParameters &EPGN);

#endif
