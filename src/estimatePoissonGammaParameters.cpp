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
EstimatePoissonGammaAlleleParameters::EstimatePoissonGammaAlleleParameters(const Eigen::VectorXd & coverage, const std::vector<Eigen::MatrixXd> & expectedContributionMatrix,
                                                                           const std::vector<Eigen::VectorXd> & alleleIndex, const Eigen::VectorXd & markerImbalances,
                                                                           const Eigen::VectorXd & partialSumAlleles, const double & convexMarkerImbalanceInterpolation,
                                                                           const Eigen::VectorXd & tolerance)
{
    Coverage = coverage;
    ExpectedContributionMatrix = expectedContributionMatrix;
    AlleleIndex = alleleIndex;

    PartialSumAlleles = partialSumAlleles;
    MarkerImbalances = markerImbalances;

    ConvexMarkerImbalanceInterpolation = convexMarkerImbalanceInterpolation;
    Tolerance = tolerance;
    MaximumNumberOfIterations = 2000;

    Counter = 0;

    NumberOfContributors = expectedContributionMatrix[0].cols();
    NumberOfMarkers = MarkerImbalances.size();

    initialiseParameters();
}


void EstimatePoissonGammaAlleleParameters::initialiseParameters()
{
    // Initialising mixture parameters
    double sumCoverage = Coverage.sum();

    Eigen::VectorXd mixtureParameters = Eigen::VectorXd::Ones(NumberOfContributors) / NumberOfContributors;
    for (std::size_t m = 0; m < NumberOfMarkers; m++)
    {
        const Eigen::VectorXd & AlleleIndex_m = AlleleIndex[m];
        const Eigen::MatrixXd & ExpectedContributionMatrix_m = ExpectedContributionMatrix[m];
        for (std::size_t a = 0; a < AlleleIndex_m.size(); a++)
        {
            std::size_t n = PartialSumAlleles[m] + AlleleIndex_m[a];
            for (std::size_t c = 0; c < NumberOfContributors; c++)
            {
                Eigen::VectorXd contributor_c = ExpectedContributionMatrix_m.col(c);
                if (contributor_c[a] >= 1)
                {
                    mixtureParameters[c] += Coverage[n] / contributor_c[a];
                }
            }
        }
    }

    MixtureParameters = mixtureParameters / mixtureParameters.sum();

    // Initialising average heterozygote parameter
    Eigen::VectorXd markerCoverage = Eigen::VectorXd::Zero(NumberOfMarkers);
    for (std::size_t m = 0; m < NumberOfMarkers; m++)
    {
        const Eigen::VectorXd & AlleleIndex_m = AlleleIndex[m];
        for (std::size_t a = 0; a < AlleleIndex_m.size(); a++)
        {
            std::size_t n = PartialSumAlleles[m] + AlleleIndex_m[a];
            markerCoverage[m] += Coverage[n] / (2 * NumberOfContributors);
        }
    }

    SampleParameters = 2.0 * Eigen::VectorXd::Ones(2);
    SampleParameters[0] = markerCoverage.sum() / markerCoverage.size();

    // Creating moment estimates of the marker imbalances
    Eigen::VectorXd markerImbalancesMoM = markerCoverage / SampleParameters[0];
    markerImbalancesMoM = markerImbalancesMoM / markerImbalancesMoM.mean();

    MarkerImbalancesParameters = ConvexMarkerImbalanceInterpolation * markerImbalancesMoM + (1 - ConvexMarkerImbalanceInterpolation) * MarkerImbalances;
}

double logLikelihoodAlleleCoverage(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    // Unpacking data
    EstimatePoissonGammaAlleleParameters *EPGA = reinterpret_cast<EstimatePoissonGammaAlleleParameters*>(data);
    const std::vector<Eigen::MatrixXd> & ExpectedContributionMatrix = EPGA->ExpectedContributionMatrix;
    const std::vector<Eigen::VectorXd> & AlleleIndex = EPGA->AlleleIndex;

    const Eigen::VectorXd & MarkerImbalances = EPGA->MarkerImbalancesParameters;
    const Eigen::VectorXd & PartialSumAlleles = EPGA->PartialSumAlleles;
    const Eigen::VectorXd & Coverage = EPGA->Coverage;

    const std::size_t NumberOfMarkers = EPGA->NumberOfMarkers;

    std::size_t & Counter = EPGA->Counter;

    const std::size_t & C = EPGA->NumberOfContributors;
    const std::size_t S = 2;

    const double & referenceMarkerAverage = x[0];
    const double & dispersion = x[1];

    double logLikelihood = 0.0;
    std::vector<double> gradient(C + S - 1, 0.0);
    for (std::size_t m = 0; m < NumberOfMarkers; m++)
    {
        const Eigen::VectorXd & AlleleIndex_m = AlleleIndex[m];
        const Eigen::MatrixXd & ExpectedContributionMatrix_m = ExpectedContributionMatrix[m];
        for (std::size_t a = 0; a < AlleleIndex_m.size(); a++)
        {
            std::size_t n = PartialSumAlleles[m] + AlleleIndex_m[a];
            const Eigen::VectorXd ExpectedContributionMatrix_ma = ExpectedContributionMatrix_m.row(a);

            double EC_n = ExpectedContributionMatrix_ma[C - 1];
            for (std::size_t c = 0; c < C - 1; c++)
            {
                EC_n += (ExpectedContributionMatrix_ma[c] - ExpectedContributionMatrix_ma[C - 1]) * x[S + c];
            }

            double mu_ma = referenceMarkerAverage * MarkerImbalances[m] * EC_n;
            if (!(mu_ma > 0.0))
            {
                mu_ma += 2e-16;
            }

            logLikelihood += logPoissonGammaDistribution(Coverage[n], mu_ma, mu_ma / dispersion);

            if (!grad.empty())
            {
                const double & gradient_common_term = boost::math::digamma(Coverage[n] + mu_ma / dispersion) - boost::math::digamma(mu_ma / dispersion) - std::log(dispersion + 1.0);

                gradient[0] += gradient_common_term * (MarkerImbalances[m] * EC_n) / (dispersion);
                gradient[1] += (Coverage[n] - mu_ma) / (dispersion * (dispersion + 1.0)) - gradient_common_term * (mu_ma / (std::pow(dispersion, 2.0)));

                for (std::size_t c = 0; c < C - 1; c++)
                {
                    const double & EC_nm = ExpectedContributionMatrix_ma[c] - ExpectedContributionMatrix_ma[C - 1];
                    gradient[S + c] +=  gradient_common_term * ((referenceMarkerAverage * MarkerImbalances[m] * EC_nm) / dispersion);
                }

            }
        }
    }

    if (!grad.empty())
    {
        grad = gradient;
    }

    Counter++;
    return logLikelihood;
}

void estimateParametersAlleleCoverage(EstimatePoissonGammaAlleleParameters &EPGA)
{
    std::size_t S = EPGA.SampleParameters.size();
    std::size_t C = EPGA.MixtureParameters.size();

    std::vector<double> parameters(S + C - 1, 0.0);
    for (std::size_t s = 0; s < S; s++)
    {
        parameters[s] = EPGA.SampleParameters[s];
    }

    for (std::size_t c = 0; c < C - 1; c++)
    {
        parameters[S + c] = EPGA.MixtureParameters[c];
    }

    // Optimiser
    // nlopt::opt individualOptimisation(nlopt::LN_BOBYQA, S + M - 1);
    nlopt::opt individualOptimisation(nlopt::LD_LBFGS, S + C - 1);

    // Box-constraints
    std::size_t N = EPGA.Coverage.size();
    std::vector<double> lowerBound(S + C - 1), upperBound(S + C - 1);
    lowerBound[0] = 1;
    lowerBound[1] = 2e-8;
    upperBound[0] = EPGA.Coverage.maxCoeff() + 1;
    upperBound[1] = std::exp(std::log(2.0) + std::log(N) - std::log(2.0 * N - 2.0))  * (std::pow(EPGA.Coverage.maxCoeff(), 2.0) + 1);

    if (C == 1)
    {
        lowerBound[S] = 1 - 2e-16;
        upperBound[S] = 1;
    }
    else
    {
        for (std::size_t j = 0; j < C - 1; j++)
        {
            lowerBound[S + j] = 2e-16;
            upperBound[S + j] = 1.0 - 2e-16;
        }
    }

    individualOptimisation.set_lower_bounds(lowerBound);
    individualOptimisation.set_upper_bounds(upperBound);

    // Objective function
    individualOptimisation.set_max_objective(logLikelihoodAlleleCoverage, &EPGA);

    individualOptimisation.set_ftol_rel(EPGA.Tolerance[0]);
    individualOptimisation.set_ftol_abs(EPGA.Tolerance[1]);
    individualOptimisation.set_xtol_rel(EPGA.Tolerance[2]);
    individualOptimisation.set_xtol_abs(EPGA.Tolerance[3]);

    individualOptimisation.set_maxeval(EPGA.MaximumNumberOfIterations);

    double logLikelihood;

    nlopt::result result = individualOptimisation.optimize(parameters, logLikelihood);

    EPGA.LogLikelihood = logLikelihood;
    Eigen::VectorXd Parameters = STDEigen(parameters);

    EPGA.SampleParameters = Parameters.segment(0, S);

    EPGA.MixtureParameters.segment(0, C - 1) = Parameters.segment(S, C - 1);
    EPGA.MixtureParameters[C - 1] =  1.0 - Parameters.segment(S, C - 1).sum();
}


// Noise coverage
EstimatePoissonGammaNoiseParameters::EstimatePoissonGammaNoiseParameters(const Eigen::VectorXd & coverage, const std::vector<Eigen::VectorXd> & noiseIndex,
                                                                         const Eigen::VectorXd & partialSumAlleles, const Eigen::VectorXd & tolerance)
{
    Coverage = coverage;
    NoiseIndex = noiseIndex;

    PartialSumAlleles = partialSumAlleles;

    NumberOfMarkers = NoiseIndex.size();
    Tolerance = tolerance;
    MaximumNumberOfIterations = 2000;

    Counter = 0;

    initialiseParameters();
}

void EstimatePoissonGammaNoiseParameters::initialiseParameters()
{
    // Eigen::VectorXd parameters = Eigen::VectorXd::Ones(2);
    //
    // double noiseCoverageSum = 0.0;
    // double noiseCoverageSize = 0.0;
    // for (std::size_t m = 0; m < NumberOfMarkers; m++)
    // {
    //     const Eigen::VectorXd & NoiseIndex_m = NoiseIndex[m];
    //     for (std::size_t a = 0; a < NoiseIndex_m.size(); a++)
    //     {
    //         std::size_t n = PartialSumAlleles[m] + NoiseIndex_m[a];
    //         noiseCoverageSum += Coverage[n];
    //         noiseCoverageSize++;
    //     }
    // }
    //
    // double averageNoiseCoverage = noiseCoverageSum / noiseCoverageSize;
    //
    // parameters[0] = averageNoiseCoverage;

    NoiseParameters = Eigen::VectorXd::Ones(2);
}

double logLikelihoodNoiseCoverage(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    EstimatePoissonGammaNoiseParameters *EPGN = reinterpret_cast<EstimatePoissonGammaNoiseParameters*>(data);
    const Eigen::VectorXd & Coverage = EPGN->Coverage;
    const std::vector<Eigen::VectorXd> NoiseIndex = EPGN->NoiseIndex;
    const Eigen::VectorXd & PartialSumAlleles = EPGN->PartialSumAlleles;
    const std::size_t & NumberOfMarkers = EPGN->NumberOfMarkers;
    std::size_t & Counter = EPGN->Counter;

    const double & mu_ma = x[0];
    const double & dispersion = x[1];
    double logeta = std::log(dispersion) - std::log(mu_ma + dispersion);

    double logLikelihood = 0.0;
    for (std::size_t m = 0; m < NumberOfMarkers; m++)
    {
        const Eigen::VectorXd & NoiseIndex_m = NoiseIndex[m];
        for (std::size_t a = 0; a < NoiseIndex_m.size(); a++)
        {
            std::size_t n = PartialSumAlleles[m] + NoiseIndex_m[a];
            logLikelihood += logPoissonGammaDistribution(Coverage[n], mu_ma, dispersion) - std::log(1 - std::exp(dispersion * logeta));
        }
    }

    Counter++;
    return logLikelihood;
}

void estimateParametersNoiseCoverage(EstimatePoissonGammaNoiseParameters &EPGN)
{
    const std::size_t & N = EPGN.NoiseParameters.size();
    std::vector<double> parameters = EigenSTD(EPGN.NoiseParameters);

    // Optimiser
    nlopt::opt individualOptimisation(nlopt::LN_BOBYQA, N);

    // Box-constraints
    std::vector<double> lowerBound(N, 2e-8), upperBound(N, HUGE_VAL);
    // lowerBound[0] = 1.0;
    // lowerBound[1] = 2e-8;

    individualOptimisation.set_lower_bounds(lowerBound);
    individualOptimisation.set_upper_bounds(upperBound);

    // Objective function
    individualOptimisation.set_max_objective(logLikelihoodNoiseCoverage, &EPGN);

    individualOptimisation.set_ftol_rel(EPGN.Tolerance[0]);
    individualOptimisation.set_ftol_abs(EPGN.Tolerance[1]);
    individualOptimisation.set_xtol_rel(EPGN.Tolerance[2]);
    individualOptimisation.set_xtol_abs(EPGN.Tolerance[3]);

    individualOptimisation.set_maxeval(EPGN.MaximumNumberOfIterations);

    double logLikelihood;
    nlopt::result result = individualOptimisation.optimize(parameters, logLikelihood);

    EPGN.LogLikelihood = logLikelihood;
    Eigen::VectorXd Parameters = STDEigen(parameters);
    EPGN.NoiseParameters = Parameters;
}
