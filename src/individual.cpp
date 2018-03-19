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
#include <boost/math/special_functions/factorials.hpp>

Individual::Individual()
{
    Fitness = -HUGE_VAL;
}

Individual::Individual(const ExperimentalSetup & ES)
{
    if ((ES.NumberOfContributors - ES.NumberOfKnownContributors) == 0)
    {
        Eigen::MatrixXd decodedProfile = decoding(EncodedProfile, ES.NumberOfAlleles, ES.NumberOfMarkers, ES.NumberOfContributors - ES.NumberOfKnownContributors);
        Eigen::MatrixXd expectedContributionProfile = GenerateExpectedContributionProfile(ES, decodedProfile);
        Eigen::VectorXd noiseProfile = GenerateNoiseProfile(ES, expectedContributionProfile);

        EstimateParameters(ES, decodedProfile, expectedContributionProfile, noiseProfile);
        CalculateFitness(ES, decodedProfile, expectedContributionProfile, noiseProfile);
    }
    else
    {
        Rcpp::stop("The only use of this contructor is when the number of contributors equals the number of known contributors under the proposed hypothesis.");
    }
}

Individual::Individual(const Eigen::VectorXd & encodedGenotype, const ExperimentalSetup & ES)
{
    EncodedProfile = encodedGenotype;

    Eigen::MatrixXd decodedProfile = decoding(EncodedProfile, ES.NumberOfAlleles, ES.NumberOfMarkers, ES.NumberOfContributors - ES.NumberOfKnownContributors);
    Eigen::MatrixXd expectedContributionProfile = GenerateExpectedContributionProfile(ES, decodedProfile);
    Eigen::VectorXd noiseProfile = GenerateNoiseProfile(ES, expectedContributionProfile);

    EstimateParameters(ES, decodedProfile, expectedContributionProfile, noiseProfile);
    CalculateFitness(ES, decodedProfile, expectedContributionProfile, noiseProfile);
}

Individual::Individual(const Eigen::VectorXd & encodedGenotype, const Eigen::VectorXd & sampleParameters, const Eigen::VectorXd & noiseParameters,
                       const Eigen::VectorXd & mixtureParameters, const Eigen::VectorXd & markerParameters, const ExperimentalSetup & ES)
{
    EncodedProfile = encodedGenotype;

    SampleParameters = sampleParameters;
    MixtureParameters = mixtureParameters;
    MarkerImbalanceParameters = markerParameters;
    NoiseParameters = noiseParameters;

    Eigen::MatrixXd decodedProfile = decoding(EncodedProfile, ES.NumberOfAlleles, ES.NumberOfMarkers, ES.NumberOfContributors - ES.NumberOfKnownContributors);
    Eigen::MatrixXd expectedContributionProfile = GenerateExpectedContributionProfile(ES, decodedProfile);
    Eigen::VectorXd noiseProfile = GenerateNoiseProfile(ES, expectedContributionProfile);

    CalculateFitness(ES, decodedProfile, expectedContributionProfile, noiseProfile);
}

Individual::Individual(const Eigen::VectorXd & encodedGenotype, const Eigen::VectorXd & sampleParameters, const Eigen::VectorXd & noiseParameters,
                       const Eigen::VectorXd & mixtureParameters, const Eigen::VectorXd & markerParameters, const double & fitness)
{
    EncodedProfile = encodedGenotype;

    SampleParameters = sampleParameters;
    MixtureParameters = mixtureParameters;
    MarkerImbalanceParameters = markerParameters;
    NoiseParameters = noiseParameters;

    Fitness = fitness;
}


double ParentStutterContribution(const std::size_t & currentAllele, const std::size_t & stutterRecursion, const double & levelOfRecursion,
                                 const Eigen::VectorXd & decodedProfile_mu, const std::vector<Eigen::MatrixXd> & potentialParents_m, const double & numberOfAlleles_m)
{
    if (stutterRecursion == 0)
    {
        return 0.0;
    }

    Eigen::MatrixXd potentialParents_ma = potentialParents_m[currentAllele];
    Eigen::VectorXd potentialParentIndex = potentialParents_ma.col(1);
    Eigen::VectorXd stutterContribution = potentialParents_ma.col(2);

    Eigen::VectorXd potentialParentContribution = Eigen::VectorXd::Zero(potentialParentIndex.size());
    if (potentialParentIndex[0] != -1)
    {
        for (std::size_t i = 0; i < potentialParentIndex.size(); i++)
        {
            double parentContribution = ParentStutterContribution(potentialParentIndex[i] - 1, stutterRecursion - 1, levelOfRecursion + 1, decodedProfile_mu, potentialParents_m, numberOfAlleles_m);
            potentialParentContribution[i] = stutterContribution[i] * (decodedProfile_mu[potentialParentIndex[i] - 1] + parentContribution) / levelOfRecursion; // ;
        }
    }

    double totalParentContribition = potentialParentContribution.sum();
    return totalParentContribition;
}


Eigen::MatrixXd Individual::GenerateExpectedContributionProfile(const ExperimentalSetup & ES, const Eigen::MatrixXd & decodedProfile)
{
    std::size_t numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;
    Eigen::VectorXd partialSumAlleles = partialSumEigen(ES.NumberOfAlleles);

    Eigen::MatrixXd profile(ES.Coverage.size(), ES.NumberOfContributors);
    if (ES.NumberOfKnownContributors > 0)
    {
        profile = bindColumns(ES.KnownProfiles, decodedProfile);
    }
    else
    {
        profile = decodedProfile;
    }

    Eigen::MatrixXd expectedContributionProfile = Eigen::MatrixXd::Zero(decodedProfile.rows(), ES.NumberOfContributors);
    for (std::size_t m = 0; m < ES.NumberOfMarkers; m++)
    {
        std::vector<Eigen::MatrixXd> potentialParents_m = ES.PotentialParents[m];
        Eigen::MatrixXd expectedContributionProfile_m = Eigen::MatrixXd::Zero(ES.NumberOfAlleles[m], ES.NumberOfContributors);
        Eigen::MatrixXd decodedProfile_m = profile.block(partialSumAlleles[m], 0, ES.NumberOfAlleles[m], ES.NumberOfContributors);
        for (std::size_t u = 0; u < ES.NumberOfContributors; u++)
        {
            Eigen::VectorXd decodedProfile_mu = decodedProfile_m.col(u);
            Eigen::VectorXd stutterContribution = Eigen::VectorXd::Zero(ES.NumberOfAlleles[m]);
            for (std::size_t a = 0; a < ES.NumberOfAlleles[m]; a++)
            {
                stutterContribution[a] = ParentStutterContribution(a, ES.LevelsOfStutterRecursion, 1, decodedProfile_mu, potentialParents_m, ES.NumberOfAlleles[m]);
            }

            expectedContributionProfile_m.col(u) = decodedProfile_mu + stutterContribution;
        }

        expectedContributionProfile.block(partialSumAlleles[m], 0, ES.NumberOfAlleles[m], ES.NumberOfContributors) = expectedContributionProfile_m;
    }

    return expectedContributionProfile;
}


Eigen::VectorXd Individual::GenerateNoiseProfile(const ExperimentalSetup & ES, const Eigen::MatrixXd & expectedContributionProfile)
{
    Eigen::VectorXd Ones = Eigen::VectorXd::Ones(expectedContributionProfile.cols());
    Eigen::VectorXd ExpectedContributionProfileRowSum = expectedContributionProfile * Ones;

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

void Individual::EstimateParameters(const ExperimentalSetup & ES, Eigen::MatrixXd decodedProfile, Eigen::MatrixXd expectedContributionProfile,
                                    Eigen::VectorXd noiseProfile)
{
    // Set-up E_c and N matrices
    Eigen::MatrixXd genoypeMatrix = bindColumns(ES.KnownProfiles, decodedProfile);
    Eigen::VectorXd alleleCount = genoypeMatrix.rowwise().sum();

    double noiseProfileSum = noiseProfile.sum();
    Eigen::VectorXd partialSumAlleles = partialSumEigen(ES.NumberOfAlleles);

    // Reducing the size of E_c, y, and y^(c)
    Eigen::MatrixXd expectedContributionProfileReduced = Eigen::MatrixXd::Zero(ES.Coverage.size() - noiseProfileSum, ES.NumberOfContributors);

    Eigen::VectorXd alleleCountReduced = Eigen::VectorXd::Zero(ES.Coverage.size() - noiseProfileSum);
    Eigen::VectorXd coverageAllele = Eigen::VectorXd::Zero(ES.Coverage.size() - noiseProfileSum);

    Eigen::VectorXd numberOfAllelesReduced = Eigen::VectorXd::Zero(ES.NumberOfAlleles.size());
    Eigen::VectorXd coverageNoise = Eigen::VectorXd::Zero(noiseProfileSum);

    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t m = 0;
    for (std::size_t n = 0; n < ES.Coverage.size(); n++)
    {
        if (partialSumAlleles[m + 1] <= n)
        {
            m++;
        }

        if (noiseProfile[n] == 0)
        {
            expectedContributionProfileReduced.row(i) = expectedContributionProfile.row(n);
            alleleCountReduced[i] = alleleCount[n];
            coverageAllele[i] = ES.Coverage[n];

            numberOfAllelesReduced[m] += 1;
            i++;
        }
        else
        {
            coverageNoise[j] = ES.Coverage[n];
            j++;
        }
    }

    // Estimating paramters
    std::vector<Eigen::VectorXd> EPGA = estimateParametersAlleleCoverage(coverageAllele, expectedContributionProfileReduced,
                                                                         ES.MarkerImbalances, numberOfAllelesReduced, alleleCountReduced,
                                                                         ES.ConvexMarkerImbalanceInterpolation, ES.Tolerance);

    // EstimatePoissonGammaAlleleParameters EPGA(coverageAllele, expectedContributionProfileReduced,
    //                                           ES.MarkerImbalances, numberOfAllelesReduced, alleleCountReduced,
    //                                           ES.ConvexMarkerImbalanceInterpolation, ES.Tolerance);

    Eigen::VectorXd EPGN = estimateParametersNoiseCoverage(coverageNoise);

    // EstimatePoissonGammaNoiseParameters EPGN(coverageNoise, ES.Tolerance);

    // estimateParametersAlleleCoverage(EPGA);
    // estimateParametersNoiseCoverage(EPGN);

    SampleParameters = EPGA[0]; // EPGA.SampleParameters;
    MixtureParameters = EPGA[1]; // EPGA.MixtureParameters;
    MarkerImbalanceParameters = EPGA[2]; // EPGA.MarkerImbalancesParameters;

    NoiseParameters = EPGN; // EPGN.NoiseParameters; //
}

void Individual::CalculateFitness(const ExperimentalSetup & ES, Eigen::MatrixXd decodedProfile, Eigen::MatrixXd expectedContributionProfile,
                                  Eigen::VectorXd noiseProfile)
{
    LogLikelihoodAlleleMarker = logLikelihoodAlleleCoverage(ES.Coverage, expectedContributionProfile, SampleParameters, MixtureParameters,
                                                            ES.NumberOfAlleles, MarkerImbalanceParameters);
    LogLikelihoodAllele = LogLikelihoodAlleleMarker.sum();

    LogLikelihoodNoiseMarker = logLikelihoodNoiseCoverage(ES.Coverage, noiseProfile, NoiseParameters[0], NoiseParameters[1], ES.NumberOfAlleles);
    LogLikelihoodNoise = LogLikelihoodNoiseMarker.sum();

    if ((ES.Theta < 0.0) | (ES.AlleleFrequencies.sum() == 0))
    {
        LogPriorGenotypeProbabilityMarker = Eigen::VectorXd::Zero(ES.NumberOfMarkers);
        LogPriorGenotypeProbability = 0.0;
    }
    else
    {
        LogPriorGenotypeProbabilityMarker = logPriorGenotypeProbability(ES.AlleleFrequencies, ES.Theta, decodedProfile, ES.AllKnownProfiles, ES.NumberOfMarkers, ES.NumberOfAlleles);
        LogPriorGenotypeProbability = LogPriorGenotypeProbabilityMarker.sum();
    }

    Fitness = LogLikelihoodAllele + LogLikelihoodNoise + LogPriorGenotypeProbability;
}

Eigen::VectorXd Individual::CalculateResiduals(const ExperimentalSetup & ES, const Eigen::MatrixXd & expectedContributionProfile,
                                               const Eigen::VectorXd & noiseProfile)
{
    const double & referenceMarkerAverage = SampleParameters[0];
    double dispersion;
    const Eigen::VectorXd & Coverage = ES.Coverage;
    const Eigen::VectorXd & MarkerImbalance = ES.MarkerImbalances;
    const Eigen::VectorXd EC = expectedContributionProfile * MixtureParameters;

    std::size_t N = Coverage.size();
    Eigen::VectorXd devianceResiduals = Eigen::VectorXd::Zero(N);
    for (std::size_t n = 0; n < N; n++)
    {
        double mu_ma;
        if (noiseProfile[n] == 0)
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

        devianceResiduals[n] = boost::math::sign(Coverage[n] - mu_ma) * std::pow(deviance_ma, 0.5);
    }

    return devianceResiduals;
}

Rcpp::List Individual::ReturnRcppList(const ExperimentalSetup & ES)
{

    Eigen::MatrixXd decodedProfile = decoding(EncodedProfile, ES.NumberOfAlleles, ES.NumberOfMarkers, ES.NumberOfContributors - ES.NumberOfKnownContributors);
    Eigen::MatrixXd expectedContributionProfile = GenerateExpectedContributionProfile(ES, decodedProfile);

    Eigen::VectorXd noiseProfile = GenerateNoiseProfile(ES, expectedContributionProfile);

    CalculateFitness(ES, decodedProfile, expectedContributionProfile, noiseProfile);

    return Rcpp::List::create(Rcpp::Named("EncodedUnknownProfiles") = EncodedProfile,
                              Rcpp::Named("DecodedUnknownProfiles") = decodedProfile,
                              Rcpp::Named("ExpectedContributionMatrix") = expectedContributionProfile,
                              Rcpp::Named("NoiseVector") = noiseProfile,
                              Rcpp::Named("Parameters") = Rcpp::List::create(
                                  Rcpp::Named("SampleParameters") = SampleParameters,
                                  Rcpp::Named("MixtureParameters") = MixtureParameters,
                                  Rcpp::Named("MarkerImbalanceParameters") = MarkerImbalanceParameters,
                                  Rcpp::Named("NoiseParameters") = NoiseParameters),
                                  Rcpp::Named("LogLikelihoodAlleleCoverage") = LogLikelihoodAlleleMarker,
                                  Rcpp::Named("LogLikelihoods") = Rcpp::NumericVector::create(LogLikelihoodAllele, LogLikelihoodNoise, LogPriorGenotypeProbability),
                                  Rcpp::Named("Fitness") = Fitness);
}

Rcpp::List Individual::ReturnRcppListSimplified()
{
    return Rcpp::List::create(Rcpp::Named("SampleParameters") = SampleParameters,
                              Rcpp::Named("MarkerImbalanceParameters") = MarkerImbalanceParameters,
                              Rcpp::Named("MixtureParameters") = MixtureParameters,
                              Rcpp::Named("NoiseParameters") = NoiseParameters,
                              Rcpp::Named("Fitness") = Fitness);
}

//[[Rcpp::export(.setupIndividual)]]
Rcpp::List setupIndividual(const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfContributors,
                           const std::size_t & numberOfKnownContributors, const Eigen::MatrixXd & knownProfiles, const Eigen::VectorXd & coverage,
                           const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents, const Eigen::VectorXd & markerImbalances,
                           const double & convexMarkerImbalanceInterpolation, const double & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                           const std::size_t & levelsOfStutterRecursion)
{

    const Eigen::MatrixXd allProfilesEmpty;
    const ExperimentalSetup ES(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allProfilesEmpty, coverage, potentialParents, markerImbalances, convexMarkerImbalanceInterpolation, tolerance, theta, alleleFrequencies, levelsOfStutterRecursion);
    Individual I(ES);

    Eigen::MatrixXd decodedProfile = decoding(I.EncodedProfile, ES.NumberOfAlleles, ES.NumberOfMarkers, ES.NumberOfContributors - ES.NumberOfKnownContributors);
    Eigen::MatrixXd expectedContributionProfile = I.GenerateExpectedContributionProfile(ES, decodedProfile);
    Eigen::VectorXd noiseProfile = I.GenerateNoiseProfile(ES, expectedContributionProfile);

    return Rcpp::List::create(Rcpp::Named("EncodedUnknownProfiles") = I.EncodedProfile,
                              Rcpp::Named("DecodedUnknownProfiles") = decodedProfile,
                              Rcpp::Named("ExpectedContributionMatrix") = expectedContributionProfile,
                              Rcpp::Named("NoiseVector") = noiseProfile,
                              Rcpp::Named("Parameters") = Rcpp::List::create(
                                  Rcpp::Named("SampleParameters") = I.SampleParameters,
                                  Rcpp::Named("MixtureParameters") = I.MixtureParameters,
                                  Rcpp::Named("MarkerImbalanceParameters") = I.MarkerImbalanceParameters,
                                  Rcpp::Named("NoiseParameters") = I.NoiseParameters),
                              Rcpp::Named("LogLikelihoodAlleleCoverage") = I.LogLikelihoodAlleleMarker,
                              Rcpp::Named("LogLikelihoods") = Rcpp::NumericVector::create(I.LogLikelihoodAllele, I.LogLikelihoodNoise, I.LogPriorGenotypeProbability),
                              Rcpp::Named("Fitness") = I.Fitness);
}
