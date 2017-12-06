#include <Rcpp.h>

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include "experimentalSetup.hpp"
#include "individual.hpp"
#include "encodingScheme.hpp"
#include "logLikelihoods.hpp"
#include "estimatePoissonGammaParameters.hpp"

#include "AuxiliaryFunctions.hpp"
#include <boost/math/special_functions/sign.hpp>

Individual::Individual()
{
    Fitness = -HUGE_VAL;
}

Individual::Individual(ExperimentalSetup ES)
{
    if ((ES.NumberOfContributors - ES.NumberOfKnownContributors) == 0)
    {
        ExpectedContributionProfile = GenerateExpectedContributionProfile(ES);
        NoiseProfile = GenerateNoiseProfile(ES);

        EstimateParameters(ES);
        CalculateResiduals(ES);
    }
    else
    {
        Rcpp::stop("The only use of this contructor is when the number of contributors equals the number of known contributors under the proposed hypothesis.");
    }
}

Individual::Individual(Eigen::VectorXd encodedGenotype, ExperimentalSetup ES)
{
    EncodedProfile = encodedGenotype;

    DecodedProfile = decoding(EncodedProfile, ES.NumberOfAlleles, ES.NumberOfMarkers, ES.NumberOfContributors - ES.NumberOfKnownContributors);

    ExpectedContributionProfile = GenerateExpectedContributionProfile(ES);

    NoiseProfile = GenerateNoiseProfile(ES);

    EstimateParameters(ES);
    CalculateResiduals(ES);
}

Individual::Individual(Eigen::VectorXd encodedGenotype, Eigen::VectorXd sampleParameters, Eigen::VectorXd noiseParameters,
                       Eigen::VectorXd mixtureParameters, ExperimentalSetup ES)
{
    EncodedProfile = encodedGenotype;

    DecodedProfile = decoding(EncodedProfile, ES.NumberOfAlleles, ES.NumberOfMarkers, ES.NumberOfContributors - ES.NumberOfKnownContributors);
    ExpectedContributionProfile = GenerateExpectedContributionProfile(ES);
    NoiseProfile = GenerateNoiseProfile(ES);

    SampleParameters = sampleParameters;
    NoiseParameters = noiseParameters;
    MixtureParameters = mixtureParameters;

    CalculateFitness(ES);
    CalculateResiduals(ES);
}

Eigen::MatrixXd Individual::GenerateExpectedContributionProfile(ExperimentalSetup ES)
{
    std::size_t numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;
    Eigen::VectorXd partialSumAlleles = partialSumEigen(ES.NumberOfAlleles);

    Eigen::MatrixXd profile(ES.Coverage.size(), ES.NumberOfContributors);
    if (ES.NumberOfKnownContributors > 0)
    {
        profile = bindColumns(ES.KnownProfiles, DecodedProfile);
    }
    else
    {
        profile = DecodedProfile;
    }

    Eigen::MatrixXd expectedContributionProfile = Eigen::MatrixXd::Zero(profile.rows(), ES.NumberOfContributors);
    for (std::size_t m = 0; m < ES.NumberOfMarkers; m++)
    {
        std::vector<Eigen::MatrixXd> potentialParents_m = ES.PotentialParents[m];
        Eigen::MatrixXd expectedContributionProfile_m = Eigen::MatrixXd::Zero(ES.NumberOfAlleles[m], ES.NumberOfContributors);
        Eigen::MatrixXd DecodedProfile_m = profile.block(partialSumAlleles[m], 0, ES.NumberOfAlleles[m], ES.NumberOfContributors);
        for (std::size_t u = 0; u < ES.NumberOfContributors; u++)
        {
            Eigen::VectorXd DecodedProfile_mu = DecodedProfile_m.col(u);
            Eigen::VectorXd potentialParentContribution_mu = Eigen::VectorXd::Zero(ES.NumberOfAlleles[m]);
            for (std::size_t a = 0; a < ES.NumberOfAlleles[m]; a++)
            {
                Eigen::MatrixXd potentialParents_ma = potentialParents_m[a];
                Eigen::VectorXd potentialParent = potentialParents_ma.col(0);
                Eigen::VectorXd stutterContribution = potentialParents_ma.col(1);

                Eigen::VectorXd potentialParentContribution_mua = Eigen::VectorXd::Zero(potentialParent.size());
                if (potentialParent.sum() != 0)
                {
                    for (std::size_t i = 0; i < potentialParent.size(); i++)
                    {
                        potentialParentContribution_mua[i] = potentialParent[i] * stutterContribution[i] * DecodedProfile_mu[i];
                    }
                }

                potentialParentContribution_mu[a] = potentialParentContribution_mua.sum();
            }

            expectedContributionProfile_m.col(u) = DecodedProfile_mu + potentialParentContribution_mu;
        }

        expectedContributionProfile.block(partialSumAlleles[m], 0, ES.NumberOfAlleles[m], ES.NumberOfContributors) = expectedContributionProfile_m;
    }

    return expectedContributionProfile;
}


Eigen::VectorXd Individual::GenerateNoiseProfile(ExperimentalSetup ES)
{
    Eigen::VectorXd Ones = Eigen::VectorXd::Ones(ExpectedContributionProfile.cols());
    Eigen::VectorXd ExpectedContributionProfileRowSum = ExpectedContributionProfile * Ones;

    std::size_t N = ExpectedContributionProfileRowSum.size();
    Eigen::VectorXd identifiedNoise = Eigen::VectorXd::Zero(N);
    for (std::size_t n = 0; n < N; n++)
    {
        if (!(ExpectedContributionProfileRowSum[n] > 0.0) & (ES.Coverage[n] > 0.0))
        {
            identifiedNoise[n] = 1.0;
        }
    }

    return identifiedNoise;
}

void Individual::EstimateParameters(ExperimentalSetup ES)
{
    EstimatePoissonGammaAlleleParameters EPGA(ES.Coverage, ExpectedContributionProfile, ES.MarkerImbalances, ES.Tolerance);
    estimateParametersAlleleCoverage(EPGA);

    SampleParameters = EPGA.SampleParameters;
    MixtureParameters = EPGA.MixtureParameters;

    LogLikelihoodAlleleMarker = logLikelihoodAlleleCoverage(ES.Coverage, ExpectedContributionProfile, SampleParameters, MixtureParameters, ES.NumberOfAlleles, ES.MarkerImbalances);
    LogLikelihoodAllele = EPGA.LogLikelihood;

    EstimatePoissonGammaNoiseParameters EPGN(ES.Coverage, NoiseProfile, ES.Tolerance);
    estimateParametersNoiseCoverage(EPGN);

    NoiseParameters = EPGN.NoiseParameters;
    LogLikelihoodNoise = EPGN.LogLikelihood;

    if ((ES.Theta < 0.0) | (ES.AlleleFrequencies.sum() == 0))
    {
        LogPriorGenotypeProbability = 0.0;
    }
    else
    {
        LogPriorGenotypeProbability = logPriorGenotypeProbability(ES.AlleleFrequencies, ES.Theta, DecodedProfile, ES.AllKnownProfiles, ES.NumberOfMarkers, ES.NumberOfAlleles);
    }

    Fitness = LogLikelihoodAllele + LogLikelihoodNoise + LogPriorGenotypeProbability;
}


void Individual::CalculateFitness(ExperimentalSetup ES)
{
    LogLikelihoodAlleleMarker = logLikelihoodAlleleCoverage(ES.Coverage, ExpectedContributionProfile, SampleParameters, MixtureParameters, ES.NumberOfAlleles, ES.MarkerImbalances);
    LogLikelihoodAllele = LogLikelihoodAlleleMarker.sum();
    LogLikelihoodNoise = logLikelihoodNoiseCoverage(ES.Coverage, NoiseProfile, NoiseParameters[0], NoiseParameters[1]);

    if ((ES.Theta < 0.0) | (ES.AlleleFrequencies.sum() == 0.0))
    {
        LogPriorGenotypeProbability = 0.0;
    }
    else
    {
        LogPriorGenotypeProbability = logPriorGenotypeProbability(ES.AlleleFrequencies, ES.Theta, DecodedProfile, ES.KnownProfiles, ES.NumberOfMarkers, ES.NumberOfAlleles);
    }

    Fitness = LogLikelihoodAllele + LogLikelihoodNoise + LogPriorGenotypeProbability;
}


void Individual::CalculateResiduals(ExperimentalSetup ES)
{
    const double & referenceMarkerAverage = SampleParameters[0];
    double dispersion;
    const Eigen::VectorXd & Coverage = ES.Coverage;
    const Eigen::VectorXd & MarkerImbalance = ES.MarkerImbalances;
    const Eigen::VectorXd EC = ExpectedContributionProfile * MixtureParameters;

    std::size_t N = Coverage.size();
    DevianceResiduals = Eigen::VectorXd::Zero(N);
    for (std::size_t n = 0; n < N; n++)
    {
        double mu_ma;
        if (NoiseProfile[n] == 0)
        {
            mu_ma = referenceMarkerAverage * MarkerImbalance[n] * EC[n];
            dispersion = SampleParameters[1];
        }
        else
        {
            mu_ma = NoiseParameters[0];
            dispersion = NoiseParameters[1];
        }

        double deviance_ma = (Coverage[n] + dispersion) * (std::log(mu_ma + dispersion) - std::log(Coverage[n] + dispersion));
        if (Coverage[n] > 0)
        {
            deviance_ma += Coverage[n] * (std::log(Coverage[n]) - std::log(mu_ma));
        }

        DevianceResiduals[n] = boost::math::sign(Coverage[n] - mu_ma) * std::pow(deviance_ma, 0.5);
    }
}

Rcpp::List Individual::ReturnRcppList()
{
    return Rcpp::List::create(Rcpp::Named("EncodedUnknownProfiles") = EncodedProfile,
                              Rcpp::Named("DecodedUnknownProfiles") = DecodedProfile,
                              Rcpp::Named("ExpectedContributionMatrix") = ExpectedContributionProfile,
                              Rcpp::Named("NoiseVector") = NoiseProfile,
                              Rcpp::Named("Parameters") = Rcpp::List::create(
                                  Rcpp::Named("SampleParameters") = SampleParameters,
                                  Rcpp::Named("MixtureParameters") = MixtureParameters,
                                  Rcpp::Named("NoiseParameters") = NoiseParameters),
                              Rcpp::Named("DevianceResiduals") = DevianceResiduals,
                              Rcpp::Named("LogLikelihoodAlleleCoverage") = LogLikelihoodAlleleMarker,
                              Rcpp::Named("LogLikelihoods") = Rcpp::NumericVector::create(LogLikelihoodAllele, LogLikelihoodNoise, LogPriorGenotypeProbability),
                              Rcpp::Named("Fitness") = Fitness);
}

Rcpp::List Individual::ReturnRcppListSimplified()
{
    return Rcpp::List::create(Rcpp::Named("SampleParameters") = SampleParameters,
                              Rcpp::Named("MixtureParameters") = MixtureParameters,
                              Rcpp::Named("Fitness") = Fitness);
}

//[[Rcpp::export(.setupIndividual)]]
Rcpp::List setupIndividual(std::size_t numberOfMarkers, Eigen::VectorXd numberOfAlleles, std::size_t numberOfContributors, std::size_t numberOfKnownContributors, Eigen::MatrixXd knownProfiles,
                           Eigen::VectorXd coverage, std::vector< std::vector < Eigen::MatrixXd > > potentialParents, Eigen::VectorXd markerImbalances,
                           double tolerance, double theta, Eigen::VectorXd alleleFrequencies) {

    Eigen::MatrixXd allProfilesEmpty;
    ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allProfilesEmpty, coverage, potentialParents, markerImbalances, tolerance, theta, alleleFrequencies);
    Individual I(ES);

    return Rcpp::List::create(Rcpp::Named("EncodedUnknownProfiles") = I.EncodedProfile,
                              Rcpp::Named("DecodedUnknownProfiles") = I.DecodedProfile,
                              Rcpp::Named("ExpectedContributionMatrix") = I.ExpectedContributionProfile,
                              Rcpp::Named("NoiseVector") = I.NoiseProfile,
                              Rcpp::Named("Parameters") = Rcpp::List::create(
                                  Rcpp::Named("SampleParameters") = I.SampleParameters,
                                  Rcpp::Named("MixtureParameters") = I.MixtureParameters,
                                  Rcpp::Named("NoiseParameters") = I.NoiseParameters),
                              Rcpp::Named("DevianceResiduals") = I.DevianceResiduals,
                              Rcpp::Named("LogLikelihoodAlleleCoverage") = I.LogLikelihoodAlleleMarker,
                              Rcpp::Named("LogLikelihoods") = Rcpp::NumericVector::create(I.LogLikelihoodAllele, I.LogLikelihoodNoise, I.LogPriorGenotypeProbability),
                              Rcpp::Named("Fitness") = I.Fitness);
}
