#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>

#include "experimentalSetup.hpp"
#include "encodingScheme.hpp"
#include "individual.hpp"
#include "population.hpp"
#include "evolutionaryAlgorithm.hpp"
#include "AuxiliaryFunctions.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp>


EvolutionaryAlgorithm::EvolutionaryAlgorithm() { }

EvolutionaryAlgorithm::EvolutionaryAlgorithm(ExperimentalSetup ES, const std::size_t populationSize, std::size_t seed)
{
    PopulationSize = populationSize;
    CurrentPopulation = InitialisePopulation(ES, seed);
}

EvolutionaryAlgorithm::EvolutionaryAlgorithm(ExperimentalSetup ES, const std::size_t populationSize,
                                             const std::size_t numberOfIterations, const std::size_t numberOfIterationsEqualMinMax, const std::size_t numberOfFittestIndividuals,
                                             const int parentSelectionWindowSize, const bool allowParentSurvival,
                                             const double crossoverProbability, const double mutationProbabilityLowerLimit, const double mutationDegreesOfFreedom,
                                             const Eigen::VectorXd mutationDecay, const double fractionFittestIndividuals,
                                             const std::size_t hillClimbingDirections, const std::size_t hillClimbingIterations,
                                             std::size_t seed)
{
    PopulationSize = populationSize;
    NumberOfFittestIndividuals = numberOfFittestIndividuals;
    NumberOfIterations = numberOfIterations;
    NumberOfIterationsEqualMinMax = numberOfIterationsEqualMinMax;

    ParentSelectionWindowSize = parentSelectionWindowSize;
    AllowParentSurvival = allowParentSurvival;

    CrossoverProbability = crossoverProbability;

    MutationProbabilityLowerLimit = mutationProbabilityLowerLimit;
    MutationDegreesOfFreedom = mutationDegreesOfFreedom;
    MutationDecay = mutationDecay;
    MutationDecay_t = MutationDecay[0];

    FittestEnsuredSurvivalFraction = fractionFittestIndividuals;

    HillClimbingDirections = hillClimbingDirections;
    HillClimbingIterations = hillClimbingIterations;

    CurrentPopulation = InitialisePopulation(ES, seed);
}

EvolutionaryAlgorithm::EvolutionaryAlgorithm(ExperimentalSetup ES, Population P,
                                             const std::size_t numberOfIterations, const std::size_t numberOfIterationsEqualMinMax, const std::size_t numberOfFittestIndividuals,
                                             const int parentSelectionWindowSize, const bool allowParentSurvival,
                                             const double crossoverProbability, const double mutationProbabilityLowerLimit, const double mutationDegreesOfFreedom,
                                             const Eigen::VectorXd mutationDecay, const double fractionFittestIndividuals,
                                             const std::size_t hillClimbingDirections, const std::size_t hillClimbingIterations)
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
    MutationDegreesOfFreedom = mutationDegreesOfFreedom;
    MutationDecay = mutationDecay;
    MutationDecay_t = MutationDecay[0];

    FittestEnsuredSurvivalFraction = fractionFittestIndividuals;

    HillClimbingDirections = hillClimbingDirections;
    HillClimbingIterations = hillClimbingIterations;
}

void EvolutionaryAlgorithm::RestructingIndividual(Individual & I, const ExperimentalSetup & ES)
{
    std::size_t NumberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;
    // Sorting unknown contributor mixtures
    Eigen::VectorXd Mixtures = I.MixtureParameters.segment(ES.NumberOfKnownContributors, NumberOfUnknownContributors);

    std::vector<int> sortedMixtures = sortedIndex(Mixtures);

    for (std::size_t n = 0; n < NumberOfUnknownContributors; n++)
    {
        I.MixtureParameters[ES.NumberOfKnownContributors + n] = Mixtures[sortedMixtures[n]];
    }

    // Sorting unknown genotypes
    Eigen::VectorXd unsortedEncodedProfile = I.EncodedProfile;
    for (std::size_t m = 0; m < ES.NumberOfMarkers; m++)
    {
        for (std::size_t c = 0; c < NumberOfUnknownContributors; c++)
        {
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

    I.DecodedProfile = decoding(I.EncodedProfile, ES.NumberOfAlleles, ES.NumberOfMarkers, NumberOfUnknownContributors);
    I.ExpectedContributionProfile = I.GenerateExpectedContributionProfile(ES);
}

Population EvolutionaryAlgorithm::InitialisePopulation(ExperimentalSetup & ES, std::size_t seed)
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

Eigen::VectorXd EvolutionaryAlgorithm::Crossover(const Individual & I, const Individual & J, std::size_t seed)
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

    return E;
}

Eigen::VectorXd EvolutionaryAlgorithm::CreateMutationProbability(const Individual & I, const ExperimentalSetup & ES)
{
    std::size_t N = I.DevianceResiduals.size();
    Eigen::VectorXd deviance = I.DevianceResiduals;
    // boost::math::normal_distribution<> devianceDistribution(0, MutationDecay_t);
    boost::math::students_t devianceDistribution(MutationDegreesOfFreedom);

    Eigen::VectorXd mutation = MutationProbabilityLowerLimit * Eigen::VectorXd::Ones(N);
    for (std::size_t n = 0; n < N; n++)
    {
        if (std::abs(deviance[n]) > MutationDecay_t)
        {
            mutation[n] = 1.0 - (1.0 - MutationProbabilityLowerLimit) * boost::math::pdf(devianceDistribution, deviance[n]) / boost::math::pdf(devianceDistribution, 0.0);
        }

        // mutation[n] = 1.0 - (1.0 - MutationProbabilityLowerLimit) * boost::math::pdf(devianceDistribution, deviance[n]) / boost::math::pdf(devianceDistribution, 0.0);
    }

    return mutation;
}

Eigen::VectorXd EvolutionaryAlgorithm::EncodeMutationProbability(Eigen::VectorXd mutationProability, Individual I, ExperimentalSetup & ES)
{
    Eigen::MatrixXd profile = I.DecodedProfile;
    std::size_t N = profile.rows();
    std::size_t numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;

    Eigen::VectorXd partialSumAlleles = partialSumEigen(ES.NumberOfAlleles);

    Eigen::VectorXd encodedMutation = Eigen::VectorXd::Zero(2 * numberOfUnknownContributors * ES.NumberOfMarkers);
    for (std::size_t i = 0; i < ES.NumberOfMarkers; i++)
    {
        Eigen::MatrixXd profile_m = profile.block(partialSumAlleles[i], 0, ES.NumberOfAlleles[i], numberOfUnknownContributors);

        Eigen::VectorXd mutation_m = mutationProability.segment(partialSumAlleles[i], ES.NumberOfAlleles[i]);

        for (std::size_t j = 0; j < numberOfUnknownContributors; j++)
        {
            std::size_t k = 2 * numberOfUnknownContributors * i + 2 * j;

            Eigen::Vector2i whichNonZero = nonZeroElementsOfMarker(profile_m.col(j));
            encodedMutation[k] = mutation_m[whichNonZero[0]];
            encodedMutation[k + 1] = mutation_m[whichNonZero[1]];
        }
    }

    return encodedMutation;
}

Individual EvolutionaryAlgorithm::Mutation(Eigen::VectorXd E, ExperimentalSetup & ES, std::size_t seed)
{
    Individual unmutatedIndividual(E, ES);
    Eigen::VectorXd I = unmutatedIndividual.EncodedProfile;

    boost::random::mt19937 rng(seed);
    boost::random::uniform_real_distribution<> uniform(0, 1);

    Eigen::VectorXd mutation = CreateMutationProbability(unmutatedIndividual, ES);

    Eigen::VectorXd encodedMutation = EncodeMutationProbability(mutation, unmutatedIndividual, ES);

    std::size_t numberOfUnknownContributors = ES.NumberOfContributors - ES.NumberOfKnownContributors;
    // Eigen::VectorXd partialSumAlleles = partialSumEigen(ES.NumberOfAlleles);
    for (std::size_t i = 0; i < ES.NumberOfMarkers; i++)
    {
        for (std::size_t j = 0; j < 2 * numberOfUnknownContributors; j++)
        {
            double mutate = uniform(rng);

            std::size_t k = 2 * numberOfUnknownContributors * i + j;
            if (mutate < encodedMutation[k])
            {
                // std::size_t n = partialSumAlleles[i] + I[k];
                // double I_mu_ma = unmutatedIndividual.SampleParameters[0] * ES.MarkerImbalances[n] *
                //     unmutatedIndividual.ExpectedContributionProfile.row(n) * unmutatedIndividual.MixtureParameters;
                //
                // std::vector<Individual> surroundings(ES.NumberOfAlleles[i] - 1);
                // Eigen::VectorXd surroundingResiduals = Eigen::VectorXd::Zero(ES.NumberOfAlleles[i] - 1);
                // for (std::size_t j = 1; j < ES.NumberOfAlleles[i]; j++)
                // {
                //     Eigen::VectorXd I_j = I;
                //
                //     // double randomDirection = uniformDirection(rng);
                //     I_j[k] = static_cast<int>(I_j[k] + j) % static_cast<int>(ES.NumberOfAlleles[i]);
                //     std::size_t m = partialSumAlleles[i] + I_j[k];
                //
                //     Individual J(I_j, unmutatedIndividual.SampleParameters, unmutatedIndividual.NoiseParameters, unmutatedIndividual.MixtureParameters, ES);
                //     double J_mu_ma = J.SampleParameters[0] * ES.MarkerImbalances[m] * J.ExpectedContributionProfile.row(m) * J.MixtureParameters;
                //
                //     RestructingIndividual(J, ES);
                //     surroundings[j - 1] = J;
                //     surroundingResiduals[j - 1] = std::abs(ES.Coverage[n] - I_mu_ma + ES.Coverage[m] - J_mu_ma);
                // }
                //
                // Eigen::MatrixXf::Index minIndex;
                // double smallestValue = surroundingResiduals.minCoeff(&minIndex);
                //
                // I = surroundings[minIndex].EncodedProfile;

                boost::random::uniform_int_distribution<> uniformMutation(1, ES.NumberOfAlleles[i] - 1);
                int mutationShift = uniformMutation(rng);

                I[k] = static_cast<int>(I[k] + mutationShift) % static_cast<int>(ES.NumberOfAlleles[i]);
            }
        }
    }

    Individual mutatedIndividual(I, ES);
    return mutatedIndividual;
}

std::size_t EvolutionaryAlgorithm::ChoosePartner(const Population & P, int currentIndividual, std::size_t seed)
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

void EvolutionaryAlgorithm::HillClimbing(Individual & I, ExperimentalSetup & ES, std::size_t seed)
{
    boost::random::mt19937 rng(seed);
    boost::random::uniform_int_distribution<> randomMarker(0, ES.NumberOfMarkers - 1);
    boost::random::uniform_int_distribution<> randomContributor(0, ES.NumberOfContributors - ES.NumberOfKnownContributors - 1);
    boost::random::uniform_int_distribution<> randomBinary(0, 1);

    Eigen::VectorXd partialSumAlleles = partialSumEigen(ES.NumberOfAlleles);
    for (std::size_t i = 0; i < HillClimbingIterations; i++)
    {
        int stepMarker = randomMarker(rng);
        int stepContributor = randomContributor(rng);
        int stepBinary = randomBinary(rng);
        std::size_t k = 2 * stepMarker * (ES.NumberOfContributors - ES.NumberOfKnownContributors) + 2 * stepContributor + stepBinary;

        // boost::random::uniform_int_distribution<> uniformDirection(1, ES.NumberOfAlleles[stepMarker] - 1);

        std::size_t n = partialSumAlleles[stepMarker] + I.EncodedProfile[k];
        double I_mu_ma = I.SampleParameters[0] * ES.MarkerImbalances[n] * I.ExpectedContributionProfile.row(n) * I.MixtureParameters;

        std::vector<Individual> surroundings(ES.NumberOfAlleles[stepMarker]);
        Eigen::VectorXd surroundingResiduals = Eigen::VectorXd::Zero(ES.NumberOfAlleles[stepMarker]);
        for (std::size_t j = 0; j < ES.NumberOfAlleles[stepMarker]; j++)
        {
            Eigen::VectorXd I_j = I.EncodedProfile;

            // double randomDirection = uniformDirection(rng);
            I_j[k] = static_cast<int>(I_j[k] + j) % static_cast<int>(ES.NumberOfAlleles[stepMarker]);
            std::size_t m = partialSumAlleles[stepMarker] + I_j[k];

            Individual J(I_j, I.SampleParameters, I.NoiseParameters, I.MixtureParameters, ES);
            double J_mu_ma = J.SampleParameters[0] * ES.MarkerImbalances[m] * J.ExpectedContributionProfile.row(m) * J.MixtureParameters;

            RestructingIndividual(J, ES);
            surroundings[j] = J;
            surroundingResiduals[j] = std::abs(ES.Coverage[n] - I_mu_ma + ES.Coverage[m] - J_mu_ma);
        }

        Eigen::MatrixXf::Index minIndex;
        double smallestValue = surroundingResiduals.minCoeff(&minIndex);

        Individual K(surroundings[minIndex].EncodedProfile, ES);
        if (K.Fitness > I.Fitness)
        {
            I = K;
        }

        // std::vector<Individual> surroundings(ES.NumberOfAlleles[stepMarker] - 1);
        // Eigen::VectorXd surroundingFitness = Eigen::VectorXd::Zero(ES.NumberOfAlleles[stepMarker] - 1);
        // for (std::size_t j = 1; j < ES.NumberOfAlleles[stepMarker]; j++)
        // {
        //     Eigen::VectorXd I_j = I.EncodedProfile;
        //
        //     // double randomDirection = uniformDirection(rng);
        //     I_j[k] = static_cast<int>(I_j[k] + j) % static_cast<int>(ES.NumberOfAlleles[stepMarker]);
        //     Individual K(I_j, ES);
        //
        //     surroundings[j - 1] = K;
        //     surroundingFitness[j - 1] = K.Fitness;
        // }
        //
        // Eigen::MatrixXf::Index maxIndex;
        // double maxValue = surroundingFitness.maxCoeff(&maxIndex);
        //
        // if (maxValue > I.Fitness)
        // {
        //     Individual K = surroundings[maxIndex];
        //     RestructingIndividual(K, ES);
        //
        //     I = K;
        // }
    }
}

Population EvolutionaryAlgorithm::SelectionCrossoverMutation(const Population & P, ExperimentalSetup & ES, std::size_t seed)
{
    boost::random::mt19937 rng(seed);
    boost::random::uniform_int_distribution<> uniformShift(0, 1e6);

    // Creating new child population
    std::vector<Individual> childPopulation(PopulationSize);
    for (std::size_t i = 0; i < PopulationSize; i++)
    {
        Individual parent = P.Individuals[i];
        // Parent partner selection
        std::size_t seedShift = uniformShift(rng);
        std::size_t partnerIndex = ChoosePartner(P, i, seedShift);

        // Crossover
        seedShift = uniformShift(rng);
        Eigen::VectorXd crossedParents = Crossover(parent, P.Individuals[partnerIndex], seedShift);

        // Mutation
        seedShift = uniformShift(rng);
        Individual mutatedChild = Mutation(crossedParents, ES, seedShift);

        if (AllowParentSurvival)
        {
            // Hill-climbing
            if ((HillClimbingDirections != 0) & (HillClimbingIterations != 0))
            {
                std::size_t seedShiftHillClimbing = uniformShift(rng);
                HillClimbing(parent, ES, seedShiftHillClimbing);
            }

            if (parent.Fitness < mutatedChild.Fitness)
            {
                // Restructuring
                RestructingIndividual(mutatedChild, ES);
                childPopulation[i] = mutatedChild;
            }
            else
            {
                childPopulation[i] = parent;
            }
        }
        else
        {
            // Hill-climbing
            if ((HillClimbingDirections != 0) & (HillClimbingIterations != 0))
            {
                std::size_t seedShiftHillClimbing = uniformShift(rng);
                HillClimbing(mutatedChild, ES, seedShiftHillClimbing);
            }

            // Restructuring
            RestructingIndividual(mutatedChild, ES);
            childPopulation[i] = mutatedChild;
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

void EvolutionaryAlgorithm::Run(ExperimentalSetup & ES, std::size_t seed, const bool & trace)
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
        for (std::size_t m = 0; m < PopulationSize; m++)
        {
            Individual I_m = C.Individuals[sortedFitness[m]];

            double F_m = I_m.Fitness;
            // bool isDuplicate = std::find(TF.begin(), TF.end(), F_m) != TF.end();
            if ((F_m > TF[0])) // & !isDuplicate)
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
                        isDuplicate = (K.cwiseAbs().sum() < 2e-16);
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
        }

        // Updates the termination counter
        std::vector<int> sortedFitness_C = sortedIndex(C.Fitness);
        if (std::abs(C.Fitness[sortedFitness_C[0]] - C.Fitness.mean()) / std::abs(C.Fitness[sortedFitness_C[0]]) < ES.Tolerance)
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
                        << "\t\t  Highest: " << C.Fitness[sortedFitness_C[0]] << "\n"
                        << "\t\t  Average: " << C.Fitness.mean() << "\n"
                        << "\t\t  Lowest: " << C.Fitness[sortedFitness_C[sortedFitness_C.size() - 1]] << "\n"
                        << "\t\tTermination counter: " << terminationCounter << " / " << NumberOfIterationsEqualMinMax << "\n";
        }
    }

    Population FittestMembers(TI);
    FittestMembersOfEntireRun = FittestMembers;
}
