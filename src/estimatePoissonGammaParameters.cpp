#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//[[Rcpp::depends(BH)]]
#include <boost/math/special_functions/digamma.hpp>

#include "experimentalSetup.hpp"
#include "logLikelihoods.hpp"
#include "individual.hpp"
#include "estimatePoissonGammaParameters.hpp"
#include "AuxiliaryFunctions.hpp"


#include "nlopt.hpp"

// Allele coverage
EstimatePoissonGammaAlleleParameters::EstimatePoissonGammaAlleleParameters(Eigen::VectorXd coverage,
                                                                           Eigen::MatrixXd expectedContributionMatrix,
                                                                           Eigen::VectorXd markerImbalances, Eigen::VectorXd numberOfAlleles,
                                                                           Eigen::VectorXd alleleCount,
                                                                           double convexMarkerImbalanceInterpolation, double tolerance)
{
    Coverage = coverage;
    ExpectedContributionMatrix = expectedContributionMatrix;
    MarkerImbalances = markerImbalances;
    NumberOfAlleles = numberOfAlleles;
    AlleleCount = alleleCount;

    ConvexMarkerImbalanceInterpolation = convexMarkerImbalanceInterpolation;

    ToleranceFRelative = tolerance;
    ToleranceEqualityConstraints = tolerance;

    MaximumNumberOfIterations = 2000;

    Counter = 0;

    NumberOfContributors = expectedContributionMatrix.cols();

    initialiseParameters();
}


void EstimatePoissonGammaAlleleParameters::initialiseParameters()
{
    // Initialising mixture parameters
    double sumCoverage = Coverage.sum();

    Eigen::VectorXd mixtureParameters = Eigen::VectorXd::Ones(ExpectedContributionMatrix.cols()) / ExpectedContributionMatrix.cols();
    for (std::size_t c = 0; c < mixtureParameters.size(); c++)
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

    Eigen::VectorXd partialSumAlleles = partialSumEigen(NumberOfAlleles);

    Eigen::VectorXd markerCoverage = Eigen::VectorXd::Zero(NumberOfAlleles.size());
    for (std::size_t m = 0; m < NumberOfAlleles.size(); m++)
    {
        for (std::size_t a = 0; a < NumberOfAlleles[m]; a++)
        {
            std::size_t n = partialSumAlleles[m] + a;
            if (AlleleCount[n] > 0) {
                markerCoverage[m] += Coverage[n] / (2 * NumberOfContributors);
            }
        }
    }

    SampleParameters = 2.0 * Eigen::VectorXd::Ones(2);
    SampleParameters[0] = markerCoverage.sum() / markerCoverage.size(); // (markerCoverage.array() / MarkerImbalances.array()).sum() / NumberOfAlleles.sum();

    // Creating moment estimates of the marker imbalances
    Eigen::VectorXd markerImbalancesMoM = markerCoverage / SampleParameters[0];
    markerImbalancesMoM = markerImbalancesMoM / markerImbalancesMoM.mean();

    MarkerImbalancesParameters = ConvexMarkerImbalanceInterpolation * markerImbalancesMoM + (1 - ConvexMarkerImbalanceInterpolation) * MarkerImbalances;
}

double logLikelihoodAlleleCoverage(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    // Unpacking data
    EstimatePoissonGammaAlleleParameters *EPGA = reinterpret_cast<EstimatePoissonGammaAlleleParameters*>(data);
    Eigen::MatrixXd & ExpectedContributionMatrix = EPGA->ExpectedContributionMatrix;
    Eigen::VectorXd & MarkerImbalances = EPGA->MarkerImbalancesParameters;
    Eigen::VectorXd & NumberOfAlleles = EPGA->NumberOfAlleles;
    Eigen::VectorXd & Coverage = EPGA->Coverage;
    std::size_t & counter = EPGA->Counter;

    std::size_t & C = EPGA->NumberOfContributors;
    std::size_t S = 2;

    double referenceMarkerAverage = x[0];
    double dispersion = x[1];

    double logLikelihood = 0.0;
    std::vector<double> gradient(C + S - 1, 0.0);
    Eigen::VectorXd partialSumAlleles = partialSumEigen(NumberOfAlleles);
    for (std::size_t m = 0; m < NumberOfAlleles.size(); m++)
    {
        for (std::size_t a = 0; a < NumberOfAlleles[m]; a++)
        {
            std::size_t n = partialSumAlleles[m] + a;
            double EC_n = ExpectedContributionMatrix.row(n)[C - 1];
            for (std::size_t c = 0; c < C - 1; c++)
            {
                EC_n += (ExpectedContributionMatrix.row(n)[c] - ExpectedContributionMatrix.row(n)[C - 1]) * x[S + c];
            }

            double mu_ma = referenceMarkerAverage * MarkerImbalances[m] * EC_n;
            logLikelihood += logPoissonGammaDistribution(Coverage[n], mu_ma, mu_ma / dispersion);

            if (!grad.empty())
            {
                double gradient_common_term = boost::math::digamma(Coverage[n] + mu_ma / dispersion) - boost::math::digamma(mu_ma / dispersion) - std::log(dispersion + 1.0);

                gradient[0] += gradient_common_term * (MarkerImbalances[m] * EC_n) / (dispersion);
                gradient[1] += (Coverage[n] - mu_ma) / (dispersion * (dispersion + 1.0)) - gradient_common_term * (mu_ma / (std::pow(dispersion, 2.0)));

                for (std::size_t c = 0; c < C - 1; c++)
                {
                    double EC_nm = ExpectedContributionMatrix.row(n)[c] - ExpectedContributionMatrix.row(n)[C - 1];
                    gradient[S + c] +=  gradient_common_term * ((referenceMarkerAverage * MarkerImbalances[m] * EC_nm) / dispersion);
                }

            }

            // if (!grad.empty())
            // {
            //     gradient[0] += Coverage[n] / referenceMarkerAverage - (Coverage[n] + dispersion) * (MarkerImbalances[m] * EC_n) / (mu_ma + dispersion);
            //     gradient[1] += boost::math::digamma(Coverage[n] + dispersion) - boost::math::digamma(dispersion) + 1.0 +
            //         std::log(dispersion) - std::log(mu_ma + dispersion) - (Coverage[n] + dispersion) / (mu_ma + dispersion);
            //
            //     for (std::size_t c = 0; c < C - 1; c++)
            //     {
            //         double EC_nm = ExpectedContributionMatrix.row(n)[c] - ExpectedContributionMatrix.row(n)[C - 1];
            //         gradient[S + c] += EC_nm * Coverage[n] / (EC_n) - (Coverage[n] + dispersion) * (referenceMarkerAverage * MarkerImbalances[m] * EC_nm) / (mu_ma + dispersion);
            //     }
            //
            // }

        }
    }

    if (!grad.empty())
    {
        grad = gradient;
    }

    counter++;
    return logLikelihood;
}

void estimateParametersAlleleCoverage(EstimatePoissonGammaAlleleParameters &EPGA)
{
    std::size_t S = EPGA.SampleParameters.size();
    std::size_t M = EPGA.MixtureParameters.size();

    std::vector<double> parameters(S + M - 1, 0.0);
    for (std::size_t s = 0; s < S; s++)
    {
        parameters[s] = EPGA.SampleParameters[s];
    }

    for (std::size_t m = 0; m < M - 1; m++)
    {
        parameters[S + m] = EPGA.MixtureParameters[m];
    }

    // Optimiser
    // nlopt::opt individualOptimisation(nlopt::LN_BOBYQA, S + M - 1);
    nlopt::opt individualOptimisation(nlopt::LD_LBFGS, S + M - 1);

    // Box-constraints
    std::size_t N = EPGA.Coverage.size();
    std::vector<double> lowerBound(S + M - 1), upperBound(S + M - 1);
    lowerBound[0] = 1;
    lowerBound[1] = 2e-8;
    upperBound[0] = EPGA.Coverage.maxCoeff() + 1;
    upperBound[1] = std::exp(std::log(2.0) + std::log(N) - std::log(2.0 * N - 2.0))  * (std::pow(EPGA.Coverage.maxCoeff(), 2.0) + 1);

    if (M == 1)
    {
        lowerBound[S] = 0;
        upperBound[S] = 1;
    }
    else
    {
        for (std::size_t j = 0; j < M - 1; j++)
        {
            lowerBound[S + j] = 2e-16;
            upperBound[S + j] = 1.0 - 2e-16;
        }
    }

    individualOptimisation.set_lower_bounds(lowerBound);
    individualOptimisation.set_upper_bounds(upperBound);

    // Objective function
    individualOptimisation.set_max_objective(logLikelihoodAlleleCoverage, &EPGA);

    individualOptimisation.set_ftol_rel(EPGA.ToleranceFRelative);
    individualOptimisation.set_maxeval(EPGA.MaximumNumberOfIterations);

    double logLikelihood;

    nlopt::result result = individualOptimisation.optimize(parameters, logLikelihood);

    EPGA.LogLikelihood = logLikelihood;
    Eigen::VectorXd Parameters = STDEigen(parameters);

    EPGA.SampleParameters = Parameters.segment(0, S);

    EPGA.MixtureParameters.segment(0, M - 1) = Parameters.segment(S, M - 1);
    EPGA.MixtureParameters[M - 1] =  1.0 - Parameters.segment(S, M - 1).sum();
}


// Noise coverage
EstimatePoissonGammaNoiseParameters::EstimatePoissonGammaNoiseParameters(Eigen::VectorXd coverage, double tolerance)
{
    Coverage = coverage;
    ToleranceFRelative = tolerance;
    MaximumNumberOfIterations = 2000;

    Counter = 0;

    initialiseParameters();
}

void EstimatePoissonGammaNoiseParameters::initialiseParameters()
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

double logLikelihoodNoiseCoverage(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    EstimatePoissonGammaNoiseParameters *EPGN = reinterpret_cast<EstimatePoissonGammaNoiseParameters*>(data);
    Eigen::VectorXd & Coverage = EPGN->Coverage;
    std::size_t & counter = EPGN->Counter;

    const double & mu_ma = x[0];
    const double & dispersion = x[1];

    double logLikelihood = 0.0;
    for (std::size_t n = 0; n < Coverage.size(); n++)
    {
        double logeta = std::log(dispersion) - std::log(mu_ma + dispersion);
        logLikelihood += logPoissonGammaDistribution(Coverage[n], mu_ma, dispersion) - std::log(1 - std::exp(dispersion * logeta));
    }

    counter++;
    return logLikelihood;
}

void estimateParametersNoiseCoverage(EstimatePoissonGammaNoiseParameters &EPGN)
{
    const std::size_t & N = EPGN.NoiseParameters.size();
    std::vector<double> parameters = EigenSTD(EPGN.NoiseParameters);

    // Optimiser
    nlopt::opt individualOptimisation(nlopt::LN_BOBYQA, N);

    // Box-constraints
    std::vector<double> lowerBound(N, 1e-16), upperBound(N, HUGE_VAL);
    lowerBound[0] = 1.0;

    individualOptimisation.set_lower_bounds(lowerBound);
    individualOptimisation.set_upper_bounds(upperBound);

    // Objective function
    individualOptimisation.set_max_objective(logLikelihoodNoiseCoverage, &EPGN);

    individualOptimisation.set_ftol_rel(EPGN.ToleranceFRelative);
    individualOptimisation.set_maxeval(EPGN.MaximumNumberOfIterations);

    double logLikelihood;
    nlopt::result result = individualOptimisation.optimize(parameters, logLikelihood);

    EPGN.LogLikelihood = logLikelihood;
    Eigen::VectorXd Parameters = STDEigen(parameters);
    EPGN.NoiseParameters = Parameters;
}
