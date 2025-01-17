// evolutionaryAlgorithm.hpp
#ifndef evolutionaryAlgorithm
#define evolutionaryAlgorithm

#include <boost/random/mersenne_twister.hpp>

#include "individual.hpp"
#include "population.hpp"
#include "experimentalSetup.hpp"

class EvolutionaryAlgorithm
{
    public:
        std::size_t PopulationSize;
        std::size_t NumberOfFittestIndividuals;
        std::size_t NumberOfIterations;
        std::size_t NumberOfIterationsEqualMinMax;
        std::size_t Iteration;

        // Parent selection
        int ParentSelectionWindowSize;
        bool AllowParentSurvival;

        // Crossover
        double CrossoverProbability;

        // Mutation
        double MutationProbabilityLowerLimit;
        std::size_t MutationIterations;

        double MutationDegreesOfFreedom;
        Eigen::VectorXd MutationDecay;
        double MutationDecay_t;

        double FittestEnsuredSurvivalFraction;
        double FitnessEnsuredSurvivalFraction;
        Eigen::VectorXd::Index FittestIndividualIndex;

        std::size_t HillClimbingIterations;

        Population CurrentPopulation;
        Population FittestMembersOfEntireRun;

        // Constructors
        EvolutionaryAlgorithm();
        EvolutionaryAlgorithm(ExperimentalSetup & ES, const std::size_t & populationSize, InitialPopulationRandomVariates & RV);

        EvolutionaryAlgorithm(ExperimentalSetup & ES, const std::size_t & populationSize,
                              const std::size_t & numberOfIterations, const std::size_t & numberOfIterationsEqualMinMax,
                              const std::size_t & numberOfFittestIndividuals,
                              const int & parentSelectionWindowSize, const bool & allowParentSurvival,
                              const double & crossoverProbability, const double & mutationProbabilityLowerLimit,
                              const std::size_t & mutationIterations, const double & mutationDegreesOfFreedom,
                              const Eigen::VectorXd & mutationDecay, const double & fractionFittestIndividuals,
                              const std::size_t & hillClimbingIterations, InitialPopulationRandomVariates & RV);

        EvolutionaryAlgorithm(ExperimentalSetup & ES, Population & P,
                              const std::size_t & numberOfIterations, const std::size_t & numberOfIterationsEqualMinMax,
                              const std::size_t & numberOfFittestIndividuals,
                              const int & parentSelectionWindowSize, const bool & allowParentSurvival,
                              const double & crossoverProbability, const double & mutationProbabilityLowerLimit,
                              const std::size_t & mutationIterations, const double & mutationDegreesOfFreedom,
                              const Eigen::VectorXd & mutationDecay, const double & fractionFittestIndividuals,
                              const std::size_t & hillClimbingIterations);

        // Functions
        void RestructingIndividual(Individual & I, const ExperimentalSetup & ES);

        Population InitialisePopulation(ExperimentalSetup & ES, InitialPopulationRandomVariates & RV);

        std::size_t ChoosePartner(const Population & P, int currentIndividual, RandomVariates & RV);

        Individual Crossover(const Individual & I, const Individual & J, const ExperimentalSetup & ES, RandomVariates & RV);

        Eigen::VectorXd CreateMutationProbability(Individual & I, const ExperimentalSetup & ES);
        void Mutation(Individual & I, const ExperimentalSetup & ES, RandomVariates & RV);

        void HillClimbing(Individual & I, ExperimentalSetup & ES, RandomVariates & RV);

        Population SelectionCrossoverMutation(const Population & P, ExperimentalSetup & ES, RandomVariates & RV);

        void Run(ExperimentalSetup & ES, RandomVariates & RV, const bool & trace);
};

#endif
