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
#include "cppoptlib/meta.h"
#include "cppoptlib/boundedproblem.h"
#include "cppoptlib/solver/lbfgsbsolver.h"

// Allele coverage
EstimatePoissonGammaAlleleParameters::EstimatePoissonGammaAlleleParameters(Eigen::VectorXd coverage, Eigen::MatrixXd expectedContributionMatrix,
                                                                           Eigen::VectorXd markerImbalances, double tolerance)
{
    Coverage = coverage;
    ExpectedContributionMatrix = expectedContributionMatrix;
    MarkerImbalances = markerImbalances;
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
	Eigen::VectorXd EC = ExpectedContributionMatrix * MixtureParameters;
	Eigen::VectorXd alleles = Eigen::VectorXd::Zero(Coverage.size());
    double averageCoverage = 0.0;
    for (std::size_t n = 0; n < Coverage.size(); n++)
    {
        if (EC[n] > 0)
        {
            averageCoverage += Coverage[n] / EC[n];
            alleles[n] = 1.0;
        }
    }

    SampleParameters = Eigen::VectorXd::Ones(2);
    SampleParameters[0] = averageCoverage / alleles.sum();
}

double logLikelihoodAlleleCoverage(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    // Unpacking data
    EstimatePoissonGammaAlleleParameters *EPGA = reinterpret_cast<EstimatePoissonGammaAlleleParameters*>(data);
    Eigen::MatrixXd & ExpectedContributionMatrix = EPGA->ExpectedContributionMatrix;
    Eigen::VectorXd & MarkerImbalances = EPGA->MarkerImbalances;
    Eigen::VectorXd & Coverage = EPGA->Coverage;
    std::size_t & counter = EPGA->Counter;

    std::size_t M = EPGA->NumberOfContributors;
    std::size_t S = 2;

    double referenceMarkerAverage = x[0];
    double dispersion = x[1];

    double logLikelihood = 0.0;
    std::vector<double> gradient(M + S - 1, 0.0);
    for (std::size_t n = 0; n < Coverage.size(); n++)
    {
        double EC_n = ExpectedContributionMatrix.row(n)[M - 1];
        for (std::size_t m = 0; m < M - 1; m++)
        {
            EC_n += (ExpectedContributionMatrix.row(n)[m] - ExpectedContributionMatrix.row(n)[M - 1]) * x[S + m];
        }

        double mu_ma = referenceMarkerAverage * MarkerImbalances[n] * EC_n;
        if (EC_n > 0.0)
        {
            logLikelihood += logPoissonGammaDistribution(Coverage[n], mu_ma, dispersion);

            if (!grad.empty())
            {
                gradient[0] += Coverage[n] / referenceMarkerAverage - (Coverage[n] + dispersion) * (MarkerImbalances[n] * EC_n) / (mu_ma + dispersion);
                gradient[1] += boost::math::digamma(Coverage[n] + dispersion) - boost::math::digamma(dispersion) + 1.0 +
                    std::log(dispersion) - std::log(mu_ma + dispersion) - (Coverage[n] + dispersion) / (mu_ma + dispersion);

                for (std::size_t m = 0; m < M - 1; m++)
                {
                    double EC_nm = ExpectedContributionMatrix.row(n)[m] - ExpectedContributionMatrix.row(n)[M - 1];
                    gradient[S + m] += EC_nm * Coverage[n] / (EC_n) - (Coverage[n] + dispersion) * (referenceMarkerAverage * MarkerImbalances[n] * EC_nm) / (mu_ma + dispersion);
                }

            }

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
    nlopt::opt individualOptimisation(nlopt::LD_LBFGS, S + M - 1);

    // Box-constraints
    std::vector<double> lowerBound(S + M - 1), upperBound(S + M - 1);
    for (std::size_t i = 0; i < S; i++)
    {
        lowerBound[i] = 1e-16;
        upperBound[i] = HUGE_VAL;
    }

    if (M == 1)
    {
        lowerBound[S] = 0;
        upperBound[S] = 1;
    }
    else
    {
        for (std::size_t j = 0; j < M - 1; j++)
        {
            lowerBound[S + j] = 1e-16;
            upperBound[S + j] = 1.0 - 1e-16;
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
EstimatePoissonGammaNoiseParameters::EstimatePoissonGammaNoiseParameters(Eigen::VectorXd coverage, Eigen::VectorXd noiseContribution, double tolerance)
{
    Coverage = coverage;
    NoiseContribution = noiseContribution;
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
        if (NoiseContribution[n] == 1)
        {
            averageNoiseCoverage += Coverage[n];
        }
    }

    parameters[0] = averageNoiseCoverage / NoiseContribution.sum();
    NoiseParameters = parameters;
}

double logLikelihoodNoiseCoverage(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    EstimatePoissonGammaNoiseParameters *EPGN = reinterpret_cast<EstimatePoissonGammaNoiseParameters*>(data);
    Eigen::VectorXd NoiseContribution = EPGN->NoiseContribution;
    Eigen::VectorXd Coverage = EPGN->Coverage;
    std::size_t & counter = EPGN->Counter;

    double mean = x[0];
    double dispersion = x[1];

    double logLikelihood = 0.0;
    std::vector<double> gradient(2, 0.0);
    for (std::size_t n = 0; n < Coverage.size(); n++)
    {
        double mu_ma = mean * NoiseContribution[n];
        if (mu_ma > 0.0)
        {
            logLikelihood += logPoissonGammaDistribution(Coverage[n], mu_ma, dispersion);

            if (!grad.empty())
            {
                gradient[0] += Coverage[n] / mu_ma - (Coverage[n] + dispersion) / (mu_ma + dispersion);
                gradient[1] += boost::math::digamma(Coverage[n] + dispersion) - boost::math::digamma(dispersion) + 1.0 +
                    std::log(dispersion) - std::log(mu_ma + dispersion) - (Coverage[n] + dispersion) / (mu_ma + dispersion);
            }
        }
    }

    if (!grad.empty())
    {
        grad = gradient;
    }

    counter++;
    return logLikelihood;
}

void estimateParametersNoiseCoverage(EstimatePoissonGammaNoiseParameters &EPGN)
{
    std::size_t N = EPGN.NoiseParameters.size();
    std::vector<double> parameters = EigenSTD(EPGN.NoiseParameters);

    // Optimiser
    nlopt::opt individualOptimisation(nlopt::LD_LBFGS, N);

    // Box-constraints
    std::vector<double> lowerBound(N), upperBound(N);
    for (std::size_t i = 0; i < N; i++)
    {
        lowerBound[i] = 1e-16;
        upperBound[i] = HUGE_VAL;
    }

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
