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

        // Reference population data
        double Theta;
        Eigen::VectorXd AlleleFrequencies;

        // Estimation input
        double Tolerance;

        ExperimentalSetup(std::size_t numberOfMarkers, Eigen::VectorXd numberOfAlleles, std::size_t numberOfContributors, std::size_t numberOfKnownContributors,
                          Eigen::MatrixXd knownProfiles, Eigen::MatrixXd allKnownProfiles,
                          Eigen::VectorXd coverage, std::vector< std::vector < Eigen::MatrixXd > > potentialParents, Eigen::VectorXd markerImbalances,
                          double tolerance, double theta, Eigen::VectorXd alleleFrequencies);

        Eigen::VectorXd GenerateUnknownGenotype(std::size_t seed);

};

#endif

