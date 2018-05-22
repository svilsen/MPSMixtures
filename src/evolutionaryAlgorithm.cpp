#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>

#include "experimentalSetup.hpp"
#include "encodingScheme.hpp"
#include "individual.hpp"
#include "population.hpp"
#include "evolutionaryAlgorithm.hpp"
#include "AuxiliaryFunctions.hpp"
#include "logLikelihoods.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp>


EvolutionaryAlgorithm::EvolutionaryAlgorithm() { }

EvolutionaryAlgorithm::EvolutionaryAlgorithm(ExperimentalSetup & ES, const std::size_t & populationSize, const std::size_t & seed)
{
    PopulationSize = populationSize;
    CurrentPopulation = InitialisePopulation(ES, seed);
}

EvolutionaryAlgorithm::EvolutionaryAlgorithm(ExperimentalSetup & ES, const std::size_t & populationSize,
                                             const std::size_t & numberOfIterations, const std::size_t & numberOfIterationsEqualMinMax,
                                             const std::size_t & numberOfFittestIndividuals,
                                             const int & parentSelectionWindowSize, const bool & allowParentSurvival,
                                             const double & crossoverProbability, const double & mutationProbabilityLowerLimit,
                                             const std::size_t & mutationIterations, const double & mutationDegreesOfFreedom,
                                             const Eigen::VectorXd & mutationDecay, const double & fractionFittestIndividuals,
                                             const std::size_t & hillClimbingIterations, const std::size_t & seed)
{
    PopulationSize = populationSize;
    NumberOfFittestIndividuals = numberOfFittestIndividuals;
    NumberOfIterations = numberOfIterations;
    NumberOfIterationsEqualMinMax = numberOfIterationsEqualMinMax;

    ParentSelectionWindowSize = parentSelectionWindowSize;
    AllowParentSurvival = allowParentSurvival;

    CrossoverProbability = crossoverProbability;

    MutationProbabilityLowerLimit = mutationProbabilityLowerLimit;
    MutationIterations = mutationIterations;

    MutationDegreesOfFreedom = mutationDegreesOfFreedom;
    MutationDecay = mutationDecay;
    MutationDecay_t = MutationDecay[0];

    FittestEnsuredSurvivalFraction = fractionFittestIndividuals;
    FitnessEnsuredSurvivalFraction = -HUGE_VAL;

    HillClimbingIterations = hillClimbingIterations;

    CurrentPopulation = InitialisePopulation(ES, seed);
    Iteration = 0;

    // Finds maximum fitness index
    Eigen::VectorXd::Index fittestIndividualIndex;
    CurrentPopulation.Fitness.maxCoeff(&fittestIndividualIndex);
    FittestIndividualIndex = fittestIndividualIndex;
}

EvolutionaryAlgorithm::EvolutionaryAlgorithm(ExperimentalSetup & ES, Population & P,
                                             const std::size_t & numberOfIterations, const std::size_t & numberOfIterationsEqualMinMax,
                                             const std::size_t & numberOfFittestIndividuals,
                                             const int & parentSelectionWindowSize, const bool & allowParentSurvival,
                                             const double & crossoverProbability, const double & mutationProbabilityLowerLimit,
                                             const std::size_t & mutationIterations, const double & mutationDegreesOfFreedom,
                                             const Eigen::VectorXd & mutationDecay, const double & fractionFittestIndividuals,
                                             const std::size_t & hillClimbingIterations)
{
    CurrentPopulation = P;

    PopulationSize = P.Individuals.size();
    NumberOfFittestIndividuals = numberOfFittestIndividuals;
    NumberOfIterations = numberOfIterations;
    NumberOfIterationsEqualMinMax = numberOfIterationsEqualMinMax;

    ParentSelectionWindowSize = parentSelectionWindowSize;
    AllowParentSurvival = allowParentSurvival;

    CrossoverProbability = crossoverProbability;

    MutationProbabilityLowerLimit = mutationProbabilityLowerLimit;
    MutationIterations = mutationIterations;
    MutationDegreesOfFreedom = mutationDegreesOfFreedom;
    MutationDecay = mutationDecay;
    MutationDecay_t = MutationDecay[0];

    FittestEnsuredSurvivalFraction = fractionFittestIndividuals;
    FitnessEnsuredSurvivalFraction = -HUGE_VAL;

    HillClimbingIterations = hillClimbingIterations;
    Iteration = 0;

    // Finds maximum fitness index
    Eigen::VectorXd::Index fittestIndividualIndex;
    CurrentPopulation.Fitness.maxCoeff(&fittestIndividualIndex);
    FittestIndividualIndex = fittestIndividualIndex;
}

void EvolutionaryAlgorithm::RestructingIndividual(Individual & I, const ExperimentalSetup & ES)
{
    const std::size_t & NumberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;
    // Sorting unknown contributor mixtures
    const Eigen::VectorXd & Mixtures = I.MixtureParameters.segment(ES.NumberOfKnownContributors, NumberOfUnknownContributors);

    const std::vector<int> & sortedMixtures = sortedIndex(Mixtures);
    for (std::size_t n = 0; n < NumberOfUnknownContributors; n++)
    {
        I.MixtureParameters[ES.NumberOfKnownContributors + n] = Mixtures[sortedMixtures[n]];
    }

    // Sorting unknown genotypes
    const Eigen::VectorXd unsortedEncodedProfile = I.EncodedProfile;
    const std::vector<Eigen::MatrixXd> unsortedReducedContributionMatrix = I.ReducedExpectedContributionMatrix;
    for (std::size_t m = 0; m < ES.NumberOfMarkers; m++)
    {
        const Eigen::MatrixXd & unsortedReducedContributionMatrix_m = unsortedReducedContributionMatrix[m];
        for (std::size_t c = 0; c < NumberOfUnknownContributors; c++)
        {
            I.ReducedExpectedContributionMatrix[m].col(c) = unsortedReducedContributionMatrix_m.col(sortedMixtures[c]);

            std::size_t k = 2 * NumberOfUnknownContributors * m + 2 * c;
            std::size_t contributorElement = 2 * NumberOfUnknownContributors * m + 2 * sortedMixtures[c];

            if (unsortedEncodedProfile[contributorElement] <= unsortedEncodedProfile[contributorElement + 1]) {
                I.EncodedProfile[k] = unsortedEncodedProfile[contributorElement];
                I.EncodedProfile[k + 1] = unsortedEncodedProfile[contributorElement + 1];
            }
            else {
                I.EncodedProfile[k] = unsortedEncodedProfile[contributorElement + 1];
                I.EncodedProfile[k + 1] = unsortedEncodedProfile[contributorElement];
            }
        }
    }
}

Population EvolutionaryAlgorithm::InitialisePopulation(ExperimentalSetup & ES, const std::size_t & seed)
{
    std::vector<Individual> I(PopulationSize);
    for (std::size_t n = 0; n < PopulationSize; n++)
    {
        Eigen::VectorXd U = ES.GenerateUnknownGenotype(seed + n);
        Individual I_n = Individual(U, ES);
        RestructingIndividual(I_n, ES);
        I[n] = I_n;
    }

    Population P(I);
    return P;
}


std::size_t EvolutionaryAlgorithm::ChoosePartner(const Population & P, int currentIndividual, const std::size_t & seed)
{
    boost::random::mt19937 rng(seed);

    int lowerWindow = std::round(std::max(0.0, static_cast<double>(currentIndividual - ParentSelectionWindowSize)));
    int upperWindow = std::round(std::min(currentIndividual + ParentSelectionWindowSize, static_cast<int>(PopulationSize - 1)));

    int leftBound = (currentIndividual - ParentSelectionWindowSize + PopulationSize) % PopulationSize;
    int rightBound = (currentIndividual + ParentSelectionWindowSize) % PopulationSize;

    Eigen::VectorXd neighbourhood = Eigen::VectorXd::Zero(2 * ParentSelectionWindowSize);
    Eigen::VectorXd logNeighbourhoodFitness = Eigen::VectorXd::Zero(2 * ParentSelectionWindowSize);
    for (std::size_t i = 0; i < 2 * ParentSelectionWindowSize; i++)
    {
        bool j = (i >= ParentSelectionWindowSize);
        int k = (leftBound + i + j) % PopulationSize;
        neighbourhood[i] = k;
        logNeighbourhoodFitness[i] = P.Fitness[k];
    }

    double logNeighbourhoodFitnessCenter = logNeighbourhoodFitness.maxCoeff();
    Eigen::VectorXd neighbourhoodFitness = Eigen::VectorXd::Zero(2 * ParentSelectionWindowSize);
    for (std::size_t i = 0; i < 2 * ParentSelectionWindowSize; i++)
    {
        neighbourhoodFitness[i] = std::exp(logNeighbourhoodFitness[i] - logNeighbourhoodFitnessCenter);
    }

    Eigen::VectorXd windowProbabilities = partialSumEigen(neighbourhoodFitness / neighbourhoodFitness.sum());

    boost::random::uniform_real_distribution<> uniform(0, 1);
    double u = uniform(rng);

    int h = 0;
    while ((u > windowProbabilities[h + 1])) {
        h++;
    }

    std::size_t partnerIndex = neighbourhood[h];
    return partnerIndex;
}

Individual EvolutionaryAlgorithm::Crossover(const Individual & I, const Individual & J, const ExperimentalSetup & ES, const std::size_t & seed)
{
    Eigen::MatrixXd E_IJ = bindColumns(I.EncodedProfile, J.EncodedProfile);

    boost::random::mt19937 rng(seed);
    boost::random::uniform_int_distribution<> uniform01(0, 1);
    boost::random::uniform_real_distribution<> uniform(0, 1);

    int columnIndex = uniform01(rng);

    std::size_t N = E_IJ.rows();
    Eigen::VectorXd E = Eigen::VectorXd::Zero(N);
    for (std::size_t n = 0; n < N; n++)
    {
        double p = uniform(rng);
        if (p < CrossoverProbability)
        {
            columnIndex = (columnIndex + 1) % 2;
        }

        E[n] = E_IJ.row(n)[columnIndex];
    }

    Individual K(E, ES);
    RestructingIndividual(K, ES);

    return K;
}

Eigen::VectorXd EvolutionaryAlgorithm::CreateMutationProbability(Individual & I, const ExperimentalSetup & ES)
{
    const Eigen::VectorXd & E = I.EncodedProfile;
    const double & referenceMarkerAverage = I.SampleParameters[0];
    const double & dispersion = I.SampleParameters[1];

    const Eigen::VectorXd & Coverage = ES.Coverage;
    const Eigen::VectorXd & MarkerImbalance = I.MarkerImbalanceParameters;
    const Eigen::VectorXd & NumberOfAlleles = ES.NumberOfAlleles;

    std::size_t N = ES.Coverage.size();
    boost::math::students_t devianceDistribution(MutationDegreesOfFreedom);

    Eigen::VectorXd mutation = MutationProbabilityLowerLimit * Eigen::VectorXd::Ones(E.size());
    if (std::abs(MutationDecay_t - MutationProbabilityLowerLimit) > DBL_EPSILON)
    {
        for (std::size_t m = 0; m < MarkerImbalance.size(); m++)
        {
            const Eigen::VectorXd & alleleIndex_m = I.ReducedAlleleIndex[m];
            const std::size_t A = alleleIndex_m.size();

            const Eigen::MatrixXd & expectedContributionProfile_m = I.ReducedExpectedContributionMatrix[m];
            const Eigen::VectorXd & EC = expectedContributionProfile_m * I.MixtureParameters;

            const Eigen::VectorXd & E_m = E.segment(2 * (ES.NumberOfContributors - ES.NumberOfKnownContributors) * m, 2 * (ES.NumberOfContributors - ES.NumberOfKnownContributors));
            for (std::size_t i = 0; i < E_m.size(); i++)
            {
                const std::size_t & n = ES.PartialSumAlleles[m] + E_m[i];
                std::size_t a = 0;
                while ((alleleIndex_m[a] != E_m[i]) & (a < A - 1))
                {
                    a++;
                }

                double mu_ma = referenceMarkerAverage * MarkerImbalance[m] * EC[a];
                if ((mu_ma == 0) & (Coverage[n] != 0))
                {
                    mu_ma += 2e-8;
                }

                double deviance_n = devianceResidualPoissonGammaDistribution(Coverage[n], mu_ma, mu_ma / dispersion);
                if (std::isnan(deviance_n))
                {
                    deviance_n = 0.0;
                    Rcpp::warning("Deviance returned 'nan'.");
                }
                else if (std::isinf(deviance_n))
                {
                    deviance_n = HUGE_VAL;
                    Rcpp::warning("Deviance returned 'inf'.");
                }

                mutation[2 * (ES.NumberOfContributors - ES.NumberOfKnownContributors) * m + i] =
                    MutationDecay_t - (MutationDecay_t - MutationProbabilityLowerLimit) * boost::math::pdf(devianceDistribution, deviance_n) /
                        boost::math::pdf(devianceDistribution, 0.0);
            }
        }
    }

    return mutation;
}

void updateMarkerIndividual(const Eigen::VectorXd & E, std::vector<Eigen::MatrixXd> reducedExpectedContributionMatrix, std::vector<Eigen::VectorXd> reducedAlleleIndex,
                            std::vector<Eigen::VectorXd> reducedNoiseIndex, const std::size_t & m, const ExperimentalSetup & ES)
{
    const std::size_t & numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;

    // Creating genotype matrix of marker m
    Eigen::MatrixXd decodedProfile_m = Eigen::MatrixXd::Zero(ES.NumberOfAlleles[m], ES.NumberOfContributors);
    for (std::size_t i = 0; i < ES.NumberOfAlleles[m]; i++)
    {
        for (std::size_t j = 0; j < ES.NumberOfKnownContributors; j++)
        {
            decodedProfile_m(i, j) = ES.KnownProfiles(ES.PartialSumAlleles[m] + i, j);
        }
    }

    for (std::size_t i = 0; i < numberOfUnknownContributors; i++)
    {
        for (std::size_t j = 0; j <= 1; j++)
        {
            std::size_t l = 2 * m * numberOfUnknownContributors + 2 * i + j;
            decodedProfile_m(E[l], i + ES.NumberOfKnownContributors) += 1;
        }
    }

    // Creating ECM of row n
    std::vector<Eigen::MatrixXd> potentialParents_m = ES.PotentialParents[m];
    Eigen::MatrixXd expectedContributionProfile_m = Eigen::MatrixXd::Zero(ES.NumberOfAlleles[m], ES.NumberOfContributors);
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

    double noiseProfileSum_m = 0.0;
    double noiseProfileSize_m = expectedContributionProfile_m.rows();
    for (std::size_t i = 0; i < noiseProfileSize_m; i++)
    {
        const double & ecps = expectedContributionProfile_m.row(i).sum();
        if (!(ecps > 2e-16))
        {
            noiseProfileSum_m++;
        }
    }

    Eigen::MatrixXd reducedExpectedContributionMatrix_m = Eigen::MatrixXd::Zero(noiseProfileSize_m - noiseProfileSum_m, ES.NumberOfContributors);
    Eigen::VectorXd reducedAlleleIndex_m = Eigen::VectorXd::Zero(noiseProfileSize_m - noiseProfileSum_m);
    Eigen::VectorXd reducedNoiseIndex_m = Eigen::VectorXd::Zero(noiseProfileSum_m);

    std::size_t n = 0;
    std::size_t i = 0, j = 0;
    for (std::size_t a = 0; a < ES.NumberOfAlleles[m]; a++)
    {
        const Eigen::VectorXd & expectedContributionProfile_mn = expectedContributionProfile_m.row(n);

        if (expectedContributionProfile_mn.sum() > 0)
        {
            reducedExpectedContributionMatrix_m.row(i) = expectedContributionProfile_mn;
            reducedAlleleIndex_m[i] = a;
            i++;
        }
        else
        {
            reducedNoiseIndex_m[j] = a;
            j++;
        }

        n++;
    }

    reducedExpectedContributionMatrix[m] = reducedExpectedContributionMatrix_m;
    reducedAlleleIndex[m] = reducedAlleleIndex_m;
    reducedNoiseIndex[m] = reducedNoiseIndex_m;
}

void EvolutionaryAlgorithm::Mutation(Individual & I, const ExperimentalSetup & ES, const std::size_t & seed)
{
    boost::random::mt19937 rng(seed);
    boost::random::uniform_real_distribution<> uniform(0, 1);

    std::size_t numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;

    Eigen::VectorXd E = I.EncodedProfile;
    Eigen::VectorXd encodedMutation = CreateMutationProbability(I, ES);

    const std::vector<Eigen::MatrixXd> reducedExpectedContributionMatrix = I.ReducedExpectedContributionMatrix;
    const std::vector<Eigen::VectorXd> reducedAlleleIndex = I.ReducedAlleleIndex;
    const std::vector<Eigen::VectorXd> reducedNoiseIndex = I.ReducedNoiseIndex;
    for (std::size_t m = 0; m < ES.NumberOfMarkers; m++)
    {
        std::size_t j = 0;
        const std::size_t i = 2 * numberOfUnknownContributors * m;
        for (std::size_t k = 0; k < 2 * numberOfUnknownContributors; k++)
        {
            double mutate = uniform(rng);
            if (mutate < encodedMutation[i + k])
            {
                boost::random::uniform_int_distribution<> uniformMutation(1, ES.NumberOfAlleles[m] - 1);
                int mutationShift = uniformMutation(rng);

                const double E_k = E[i + k];
                E[i + k] = static_cast<int>(E_k + mutationShift) % static_cast<int>(ES.NumberOfAlleles[m]);
                j++;
            }
        }

        if (j > 0)
        {
            updateMarkerIndividual(E, reducedExpectedContributionMatrix, reducedAlleleIndex, reducedNoiseIndex, m, ES);
        }
    }

    Individual J(E, reducedExpectedContributionMatrix, reducedAlleleIndex, reducedNoiseIndex, ES);
    I = J;
}

Eigen::VectorXd expectedContributionMatrixRow(const Eigen::VectorXd & E, const ExperimentalSetup & ES,
                                              const int & m, const std::size_t & k, const std::size_t & n)
{

    const std::size_t & a = E[k];
    const std::size_t & numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;

    // Creating genotype matrix of marker m
    Eigen::MatrixXd decodedProfile_m = Eigen::MatrixXd::Zero(ES.NumberOfAlleles[m], ES.NumberOfContributors);
    for (std::size_t i = 0; i < ES.NumberOfAlleles[m]; i++)
    {
        for (std::size_t j = 0; j < ES.NumberOfKnownContributors; j++)
        {
            decodedProfile_m(i, j) = ES.KnownProfiles(ES.PartialSumAlleles[m] + i, j);
        }
    }

    for (std::size_t i = 0; i < numberOfUnknownContributors; i++)
    {
        for (std::size_t j = 0; j <= 1; j++)
        {
            std::size_t l = 2 * m * numberOfUnknownContributors + 2 * i + j;
            decodedProfile_m(E[l], i + ES.NumberOfKnownContributors) += 1;
        }
    }

    // Creating ECM of row n
    std::vector<Eigen::MatrixXd> potentialParents_m = ES.PotentialParents[m];
    Eigen::VectorXd ECM_n = Eigen::VectorXd::Zero(decodedProfile_m.cols());
    for (std::size_t j = 0; j < ECM_n.size(); j++)
    {
        const Eigen::VectorXd & decodedProfile_mj = decodedProfile_m.col(j);

        double potentialParentContribution = ParentStutterContribution(a, ES.LevelsOfStutterRecursion, 1, decodedProfile_mj, potentialParents_m, ES.NumberOfAlleles[m]);
        ECM_n[j] = decodedProfile_mj[a] + potentialParentContribution;
    }

    return ECM_n;
}

void EvolutionaryAlgorithm::HillClimbing(Individual & I, ExperimentalSetup & ES, const std::size_t & seed)
{
    std::size_t numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;
    boost::random::mt19937 rng(seed);
    boost::random::uniform_int_distribution<> randomMarker(0, ES.NumberOfMarkers - 1);
    boost::random::uniform_int_distribution<> randomContributor(0, numberOfUnknownContributors - 1);
    boost::random::uniform_int_distribution<> randomBinary(0, 1);

    for (std::size_t i = 0; i < HillClimbingIterations; i++)
    {
        int stepMarker = randomMarker(rng);
        int stepContributor = randomContributor(rng);
        int stepBinary = randomBinary(rng);
        std::size_t k = 2 * stepMarker * numberOfUnknownContributors + 2 * stepContributor + stepBinary;

        // Create n'th row of the decoded and ecm matrix
        std::size_t n = ES.PartialSumAlleles[stepMarker] + I.EncodedProfile[k];
        Eigen::VectorXd expectedContributionMatrixRow_n_i = expectedContributionMatrixRow(I.EncodedProfile, ES, stepMarker, k, n);

        double I_mu_ma = I.SampleParameters[0] * ES.MarkerImbalances[n] * expectedContributionMatrixRow_n_i.transpose() * I.MixtureParameters;

        std::vector< Eigen::VectorXd > surroundings(ES.NumberOfAlleles[stepMarker] - 1);
        Eigen::VectorXd surroundingResiduals = Eigen::VectorXd::Zero(ES.NumberOfAlleles[stepMarker] - 1);
        for (std::size_t j = 1; j < ES.NumberOfAlleles[stepMarker]; j++)
        {
            Eigen::VectorXd I_j = I.EncodedProfile;

            I_j[k] = static_cast<int>(I_j[k] + j) % static_cast<int>(ES.NumberOfAlleles[stepMarker]);
            std::size_t m = ES.PartialSumAlleles[stepMarker] + I_j[k];

            Eigen::VectorXd expectedContributionMatrixRow_n_j = expectedContributionMatrixRow(I_j, ES, stepMarker, k, m);
            double J_mu_ma = I.SampleParameters[0] * ES.MarkerImbalances[m] * expectedContributionMatrixRow_n_j.transpose() * I.MixtureParameters;

            surroundings[j - 1] = I_j;
            surroundingResiduals[j - 1] = std::abs(ES.Coverage[n] - I_mu_ma + ES.Coverage[m] - J_mu_ma);
        }

        Eigen::MatrixXf::Index minIndex;
        double smallestValue = surroundingResiduals.minCoeff(&minIndex);

        std::vector<Eigen::MatrixXd> reducedExpectedContributionMatrix = I.ReducedExpectedContributionMatrix;
        std::vector<Eigen::VectorXd> reducedAlleleIndex = I.ReducedAlleleIndex;
        std::vector<Eigen::VectorXd> reducedNoiseIndex = I.ReducedNoiseIndex;

        updateMarkerIndividual(surroundings[minIndex], reducedExpectedContributionMatrix, reducedAlleleIndex, reducedNoiseIndex, stepMarker, ES);
        Individual K(surroundings[minIndex], reducedExpectedContributionMatrix, reducedAlleleIndex, reducedNoiseIndex, ES);
        if (K.Fitness > I.Fitness)
        {
            I = K;
        }
    }
}

Population EvolutionaryAlgorithm::SelectionCrossoverMutation(const Population & P, ExperimentalSetup & ES, const std::size_t & seed)
{
    boost::random::mt19937 rng(seed);
    boost::random::uniform_int_distribution<> uniformShift(0, 1e6);

    // Creating new child population
    std::vector<Individual> childPopulation(PopulationSize);
    for (std::size_t i = 0; i < PopulationSize; i++)
    {
        Individual parent = P.Individuals[i];

        // Parent partner selection
        std::size_t seedShift_1 = uniformShift(rng);
        std::size_t partnerIndex = ChoosePartner(P, i, seedShift_1);

        // Crossover
        std::size_t seedShift_2 = uniformShift(rng);
        Individual child = Crossover(parent, P.Individuals[partnerIndex], ES, seedShift_2);

        // Mutation
        std::size_t seedShift_3 = uniformShift(rng);
        Mutation(child, ES, seedShift_3);

        if (AllowParentSurvival)
        {
            // Hill-climbing
            if (HillClimbingIterations != 0)
            {
                std::size_t seedShiftHillClimbing = uniformShift(rng);
                HillClimbing(parent, ES, seedShiftHillClimbing);
                RestructingIndividual(parent, ES);
            }

            if (parent.Fitness < child.Fitness)
            {
                // Restructuring
                RestructingIndividual(child, ES);
                childPopulation[i] = child;
            }
            else
            {
                childPopulation[i] = parent;
            }
        }
        else
        {
            // Hill-climbing
            if (HillClimbingIterations != 0)
            {
                std::size_t seedShiftHillClimbing = uniformShift(rng);
                HillClimbing(child, ES, seedShiftHillClimbing);
            }

            // Restructuring
            RestructingIndividual(child, ES);
            childPopulation[i] = child;
        }


    }

    Population C(childPopulation);
    return C;
}

bool individualsEqual(Individual & I, Individual & J)
{
    Eigen::VectorXd K = I.EncodedProfile - J.EncodedProfile;
    return (K.sum() < 2e-16);
}

void EvolutionaryAlgorithm::Run(ExperimentalSetup & ES, const std::size_t & seed, const bool & trace)
{
    boost::random::mt19937 rng(seed);
    boost::random::uniform_int_distribution<> uniformShift(1, 1e6);

    std::vector<Individual> TI(NumberOfFittestIndividuals);
    std::vector<double> TF(NumberOfFittestIndividuals, -HUGE_VAL);

    std::size_t n = 0;
    std::size_t terminationCounter = 0;
    bool terminate = false;
    while (!terminate)
    {
        // Updating current population
        MutationDecay_t = MutationDecay[n];

        Population C = SelectionCrossoverMutation(CurrentPopulation, ES, seed + uniformShift(rng));
        CurrentPopulation = C;

        // Updating the list of fittest individuals
        std::vector<int> sortedFitness = sortedIndex(C.Fitness);

        std::size_t m = 0;
        bool updateFittest = true;
        while (updateFittest)
        {
            Individual I_m = C.Individuals[sortedFitness[m]];

            double F_m = I_m.Fitness;
            if (F_m > TF[0])
            {
                auto it_front = std::lower_bound(TF.cbegin(), TF.cend(), F_m);
                auto it_end = std::upper_bound(TF.cbegin(), TF.cend(), F_m);
                std::size_t i = it_front - TF.cbegin();
                std::size_t j = it_end - TF.cbegin();

                bool isDuplicate = false;
                if (F_m == TF[j - 1]) {
                    for (std::size_t k = i - 1; k < j; k++)
                    {
                        Eigen::VectorXd K = I_m.EncodedProfile - TI[k].EncodedProfile;
                        double sumAbs = 0.0;
                        for (std::size_t h = 0; h < K.size(); h++) {
                            sumAbs += std::abs(K[h]);
                        }

                        isDuplicate = (sumAbs < 2e-16);
                    }
                }

                if (!isDuplicate)
                {
                    for (std::size_t k = 0; k < j - 1; k++)
                    {
                        TI[k] = TI[k + 1];
                        TF[k] = TF[k + 1];
                    }

                    TI[j - 1] = I_m;
                    TF[j - 1] = F_m;
                }
            }

            m++;
            updateFittest = (m < NumberOfFittestIndividuals) & (F_m > TF[0]);
        }

        // Updates the termination counter
        if (std::abs(C.Fitness[sortedFitness[0]] - C.Fitness.mean()) / std::abs(C.Fitness[sortedFitness[0]]) < ES.Tolerance[0])
        {
            terminationCounter++;
        }
        else
        {
            terminationCounter = 0;
        }

        n++;
        terminate = (terminationCounter >= NumberOfIterationsEqualMinMax) || (n >= NumberOfIterations);

        if (trace)
        {
            Rcpp::Rcout << "\tCurrent iteration: " << n << "\n"
                        << "\t\tPopulation Fitness: " << "\n"
                        << "\t\t  Highest: " << C.Fitness[sortedFitness[0]] << "\n"
                        << "\t\t  Average: " << C.Fitness.mean() << "\n"
                        << "\t\t  Lowest: " << C.Fitness[sortedFitness[sortedFitness.size() - 1]] << "\n"
                        << "\t\tTermination counter: " << terminationCounter << " / " << NumberOfIterationsEqualMinMax << "\n";
        }
    }

    Population FittestMembers(TI);
    FittestMembersOfEntireRun = FittestMembers;
}
