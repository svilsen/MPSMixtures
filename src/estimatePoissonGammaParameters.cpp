#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//[[Rcpp::depends(BH)]]
#include <boost/math/special_functions/digamma.hpp>

//[[Rcpp::depends(RcppNumericalSolvers)]]
#include <cppoptlib/meta.h>
#include <cppoptlib/boundedproblem.h>
#include <cppoptlib/solver/lbfgssolver.h>
#include <cppoptlib/solver/lbfgsbsolver.h>

#include "logLikelihoods.hpp"
#include "AuxiliaryFunctions.hpp"

/////
class PoissonGammaAlleleCoverage : public cppoptlib::BoundedProblem<double> {
public :
    using Superclass = BoundedProblem<double>;
    using Vector = typename Superclass::TVector;
    using Matrix = typename Superclass::THessian;

    const Vector Coverage;
    const Vector MarkerImbalances;
    const Matrix ExpectedContributionMatrix;
    const Vector NumberOfAlleles;
    const Vector AlleleCount;

    Vector PartialSumAlleles;

    Vector SampleParameters;
    Vector MixtureParameters;
    Vector MarkerImbalancesParameters;

    const double ConvexMarkerImbalanceInterpolation;

    std::size_t Counter;

    PoissonGammaAlleleCoverage(const Eigen::VectorXd & coverage, const Eigen::MatrixXd & expectedContributionMatrix,
                               const Eigen::VectorXd & markerImbalances, const Eigen::VectorXd & numberOfAlleles,
                               const Eigen::VectorXd & alleleCount, const double & convexMarkerImbalanceInterpolation,
                               const double & tolerance) :
        Superclass(coverage.size()), Coverage(coverage), ExpectedContributionMatrix(expectedContributionMatrix),
        MarkerImbalances(markerImbalances), NumberOfAlleles(numberOfAlleles), AlleleCount(alleleCount),
        ConvexMarkerImbalanceInterpolation(convexMarkerImbalanceInterpolation),
        Counter(0) {};

    void initialiseParameters();

    double value(const Vector &theta);
    void gradient(const Vector &theta, Vector &grad);

};

void PoissonGammaAlleleCoverage::initialiseParameters()
{
    PartialSumAlleles = partialSumEigen(NumberOfAlleles);

    // Initialising mixture parameters
    double sumCoverage = Coverage.sum();
    std::size_t numberOfContributors = ExpectedContributionMatrix.cols();

    Eigen::VectorXd mixtureParameters = Eigen::VectorXd::Ones(numberOfContributors) / numberOfContributors;
    for (std::size_t c = 0; c < numberOfContributors; c++)
    {
        Eigen::VectorXd contributor_c = ExpectedContributionMatrix.col(c);

        double ECC = 0.0;
        for (std::size_t n = 0; n < Coverage.size(); n++)
        {
            if (contributor_c[n] > 0)
            {
                ECC += Coverage[n] / contributor_c[n];
            }
        }

        mixtureParameters[c] = ECC / sumCoverage;
    }
    MixtureParameters = mixtureParameters / mixtureParameters.sum();

    // Initialising average heterozygote parameter
    Eigen::VectorXd markerCoverage = Eigen::VectorXd::Zero(NumberOfAlleles.size());
    for (std::size_t m = 0; m < NumberOfAlleles.size(); m++)
    {
        for (std::size_t a = 0; a < NumberOfAlleles[m]; a++)
        {
            std::size_t n = PartialSumAlleles[m] + a;
            if (AlleleCount[n] > 0) {
                markerCoverage[m] += Coverage[n] / (2 * numberOfContributors);
            }
        }
    }

    SampleParameters = 2.0 * Eigen::VectorXd::Ones(2);
    SampleParameters[0] = markerCoverage.sum() / markerCoverage.size();

    // Creating moment estimates of the marker imbalances
    Eigen::VectorXd markerImbalancesMoM = markerCoverage / SampleParameters[0];
    markerImbalancesMoM = markerImbalancesMoM / markerImbalancesMoM.mean();

    MarkerImbalancesParameters = ConvexMarkerImbalanceInterpolation * markerImbalancesMoM + (1 - ConvexMarkerImbalanceInterpolation) * MarkerImbalances;
}

double PoissonGammaAlleleCoverage::value(const Eigen::VectorXd &theta)
{
    const std::size_t S = 2;
    const std::size_t C = theta.size() - S;

    double referenceMarkerAverage = theta[0];
    double dispersion = theta[1];

    double logLikelihood = 0.0;
    for (std::size_t m = 0; m < NumberOfAlleles.size(); m++)
    {
        for (std::size_t a = 0; a < NumberOfAlleles[m]; a++)
        {
            std::size_t n = PartialSumAlleles[m] + a;
            double EC_n = ExpectedContributionMatrix.row(n)[C];
            for (std::size_t c = 0; c < C; c++)
            {
                EC_n += (ExpectedContributionMatrix.row(n)[c] - ExpectedContributionMatrix.row(n)[C]) * theta[S + c];
            }

            double mu_ma = referenceMarkerAverage * MarkerImbalancesParameters[m] * EC_n;
            logLikelihood -= logPoissonGammaDistribution(Coverage[n], mu_ma, mu_ma / dispersion);
        }
    }

    Counter += 1;
    return logLikelihood;
}


void PoissonGammaAlleleCoverage::gradient(const Eigen::VectorXd &theta, Eigen::VectorXd &grad)
{
    std::size_t S = 2;
    std::size_t C = theta.size() - S;

    double referenceMarkerAverage = theta[0];
    double dispersion = theta[1];

    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(C + S);
    for (std::size_t m = 0; m < NumberOfAlleles.size(); m++)
    {
        for (std::size_t a = 0; a < NumberOfAlleles[m]; a++)
        {
            std::size_t n = PartialSumAlleles[m] + a;
            double EC_n = ExpectedContributionMatrix.row(n)[C];

            for (std::size_t c = 0; c < C; c++)
            {
                EC_n += (ExpectedContributionMatrix.row(n)[c] - ExpectedContributionMatrix.row(n)[C]) * theta[S + c];
            }

            double mu_ma = referenceMarkerAverage * MarkerImbalancesParameters[m] * EC_n;
            double gradient_common_term = boost::math::digamma(Coverage[n] + mu_ma / dispersion) - boost::math::digamma(mu_ma / dispersion) - std::log(dispersion + 1.0);

            gradient[0] += gradient_common_term * (MarkerImbalancesParameters[m] * EC_n) / (dispersion);
            gradient[1] += (Coverage[n] - mu_ma) / (dispersion * (dispersion + 1.0)) - gradient_common_term * (mu_ma / (std::pow(dispersion, 2.0)));

            for (std::size_t c = 0; c < C; c++)
            {
                double EC_nm = ExpectedContributionMatrix.row(n)[c] - ExpectedContributionMatrix.row(n)[C];
                gradient[S + c] +=  gradient_common_term * ((referenceMarkerAverage * MarkerImbalancesParameters[m] * EC_nm) / dispersion);
            }
        }
    }

    grad = -gradient;
}

std::vector<Eigen::VectorXd> estimateParametersAlleleCoverage(const Eigen::VectorXd & coverage, const Eigen::MatrixXd & expectedContributionMatrix,
                                                              const Eigen::VectorXd & markerImbalances, const Eigen::VectorXd & numberOfAlleles,
                                                              const Eigen::VectorXd & alleleCount, const double & convexMarkerImbalanceInterpolation,
                                                              const double & tolerance)
{
    PoissonGammaAlleleCoverage PGAC(coverage, expectedContributionMatrix, markerImbalances, numberOfAlleles,
                                    alleleCount, convexMarkerImbalanceInterpolation, tolerance);

    PGAC.initialiseParameters();

    std::size_t S = PGAC.SampleParameters.size();
    std::size_t C = PGAC.MixtureParameters.size();

    // Vector to be optimised
    Eigen::VectorXd parameters = Eigen::VectorXd::Zero(S + C - 1);
    parameters.segment(0, S) = PGAC.SampleParameters;
    parameters.segment(S, C - 1) = PGAC.MixtureParameters.segment(0, C - 1);

    // Create and set bounds
    std::size_t N = coverage.size();
    Eigen::VectorXd lowerBound(S + C - 1), upperBound(S + C - 1);
    lowerBound[0] = 1;
    lowerBound[1] = 2e-8;
    upperBound[0] = PGAC.Coverage.maxCoeff() + 1;
    upperBound[1] = std::exp(std::log(2.0) + std::log(N) - std::log(2.0 * N - 2.0)) * (std::pow(PGAC.Coverage.maxCoeff(), 2.0) + 1);

    for (std::size_t j = 0; j < C - 1; j++)
    {
        lowerBound[S + j] = 2e-8;
        upperBound[S + j] = 1.0 - 2e-8;
    }

    PGAC.setLowerBound(lowerBound);
    PGAC.setUpperBound(upperBound);

    // Setting solver and minimisation criteria
    cppoptlib::LbfgsbSolver<PoissonGammaAlleleCoverage> solver;
    cppoptlib::Criteria<double> criteria = cppoptlib::Criteria<double>::defaults();

    // criteria.fDelta = tolerance;
    // criteria.iterations = 2000;

    solver.setStopCriteria(criteria);
    solver.minimize(PGAC, parameters);

    for (std::size_t s = 0; s < S; s++)
    {
        PGAC.SampleParameters[s] = parameters[s];
    }

    double sumMixtureParameters = 0.0;
    for (std::size_t c = 0; c < C - 1; c++)
    {
        PGAC.MixtureParameters[c] = parameters[S + c];
        sumMixtureParameters += PGAC.MixtureParameters[c];
    }

    PGAC.MixtureParameters[C - 1] = 1 - sumMixtureParameters;

    // Set-up return type
    std::vector<Eigen::VectorXd> returnList(3);
    returnList[0] = PGAC.SampleParameters;
    returnList[1] = PGAC.MixtureParameters;
    returnList[2] = PGAC.MarkerImbalancesParameters;

    return returnList;
}


///////
class PoissonGammaNoiseCoverage : public cppoptlib::BoundedProblem<double> {
public :
    using Superclass = BoundedProblem<double>;
    using Vector = typename Superclass::TVector;

    const Vector Coverage;
    std::size_t Counter;

    Vector NoiseParameters;

    PoissonGammaNoiseCoverage(const Eigen::VectorXd & coverage) :
        Superclass(coverage.size()), Coverage(coverage), Counter(0) {};

    void initialiseParameters();
    double value(const Vector &theta);
};


void PoissonGammaNoiseCoverage::initialiseParameters()
{
    Eigen::VectorXd parameters = Eigen::VectorXd::Ones(2);

    double averageNoiseCoverage = 0.0;
    for (std::size_t n = 0; n < Coverage.size(); n++)
    {
        averageNoiseCoverage += Coverage[n] / Coverage.size();
    }

    parameters[0] = averageNoiseCoverage;
    NoiseParameters = parameters;
}

double PoissonGammaNoiseCoverage::value(const Eigen::VectorXd &theta)
{
    const double & mu_ma = std::exp(theta[0]) + 1.0;
    const double & dispersion = std::exp(theta[1]);

    double logLikelihood = 0.0;
    for (std::size_t n = 0; n < Coverage.size(); n++)
    {
        double logeta = std::log(dispersion) - std::log(mu_ma + dispersion);
        logLikelihood -= logPoissonGammaDistribution(Coverage[n], mu_ma, dispersion) - std::log(1 - std::exp(dispersion * logeta));
    }

    Counter++;
    return logLikelihood;
}

Eigen::VectorXd estimateParametersNoiseCoverage(const Eigen::VectorXd & coverage)
{
    PoissonGammaNoiseCoverage PGNC(coverage);
    PGNC.initialiseParameters();

    std::size_t N = PGNC.NoiseParameters.size();

    // Vector to be optimised
    Eigen::VectorXd parameters = PGNC.NoiseParameters;

    // Setting solver and minimisation criteria
    cppoptlib::LbfgsSolver<PoissonGammaNoiseCoverage> solver;
    cppoptlib::Criteria<double> criteria = cppoptlib::Criteria<double>::defaults();

    solver.setStopCriteria(criteria);
    solver.minimize(PGNC, parameters);

    Eigen::VectorXd noiseParameters(2);
    noiseParameters << (std::exp(parameters[0]) + 1), std::exp(parameters[1]);
    return noiseParameters;
}
