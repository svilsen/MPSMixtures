// experimentalSetup.hpp
#ifndef experimentalSetup
#define experimentalSetup

#include <RcppEigen.h>
#include <boost/random/mersenne_twister.hpp>

typedef Eigen::Matrix<unsigned long long, Eigen::Dynamic, 1> VectorXull;

class ExperimentalSetup
{
    public:
        // Sample data
        std::size_t NumberOfMarkers;
        VectorXull NumberOfAlleles;

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
        double ConvexMarkerImbalanceInterpolation;
        std::size_t LevelsOfStutterRecursion;


        ExperimentalSetup(const std::size_t & numberOfMarkers, const VectorXull & numberOfAlleles, const std::size_t & numberOfContributors,
                          const std::size_t & numberOfKnownContributors, const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                          const Eigen::VectorXd & coverage, const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents,
                          const Eigen::VectorXd & markerImbalances, const double & convexMarkerImbalanceInterpolation, const double & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                          const std::size_t & levelsOfStutterRecursion);

        Eigen::VectorXd GenerateUnknownGenotype(const std::size_t & seed);

};

#endif

