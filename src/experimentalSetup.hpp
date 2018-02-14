// experimentalSetup.hpp
#ifndef experimentalSetup
#define experimentalSetup

#include <RcppEigen.h>
#include <boost/random/mersenne_twister.hpp>

class ExperimentalSetup
{
    public:
        // Sample data
        std::size_t NumberOfMarkers;
        Eigen::VectorXd NumberOfAlleles;

        std::size_t NumberOfContributors;
        Eigen::MatrixXd KnownProfiles;
        std::size_t NumberOfKnownContributors;

        Eigen::MatrixXd AllKnownProfiles;

        Eigen::VectorXd Coverage;
        std::vector< std::vector< Eigen::MatrixXd > > PotentialParents;
        Eigen::VectorXd MarkerImbalances;

        std::size_t LevelsOfStutterRecursion;

        // Reference population data
        double Theta;
        Eigen::VectorXd AlleleFrequencies;

        // Estimation input
        double Tolerance;

        ExperimentalSetup(const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfContributors,
                          const std::size_t & numberOfKnownContributors, const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                          const Eigen::VectorXd & coverage, const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents,
                          const Eigen::VectorXd & markerImbalances, const double & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                          const std::size_t & levelsOfStutterRecursion);

        Eigen::VectorXd GenerateUnknownGenotype(const std::size_t & seed);

};

#endif

