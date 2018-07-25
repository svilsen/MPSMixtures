// experimentalSetup.hpp
#ifndef experimentalSetup
#define experimentalSetup

#include <RcppEigen.h>

#include <ctime>

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

struct RandomVariates
{
    typedef boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > interval_uniform_generator;

    std::size_t M;
    boost::mt19937 rng;

    boost::random::uniform_real_distribution<> uniform_real;
    boost::random::uniform_int_distribution<> uniform_01;
    boost::random::uniform_int_distribution<> uniform_marker;
    boost::random::uniform_int_distribution<> uniform_unknown_contributor;

    boost::variate_generator<boost::mt19937 &, boost::random::uniform_real_distribution<> > generate_uniform_real;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_binary;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_marker;
    boost::variate_generator<boost::mt19937 &, boost::random::uniform_int_distribution<> > generate_uniform_unknown_contributor;
    std::vector< interval_uniform_generator > generate_uniform_mutation;

    RandomVariates(const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfUnknownContributors, const std::size_t & seed);
};

class ExperimentalSetup
{
    public:
        // Sample data
        std::size_t NumberOfMarkers;
        Eigen::VectorXd NumberOfAlleles;
        Eigen::VectorXd PartialSumAlleles;

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
        Eigen::VectorXd Tolerance;
        double ConvexMarkerImbalanceInterpolation;
        std::size_t LevelsOfStutterRecursion;

        ExperimentalSetup(const std::size_t & numberOfMarkers, const Eigen::VectorXd & numberOfAlleles, const std::size_t & numberOfContributors,
                          const std::size_t & numberOfKnownContributors, const Eigen::MatrixXd & knownProfiles, const Eigen::MatrixXd & allKnownProfiles,
                          const Eigen::VectorXd & coverage, const std::vector< std::vector < Eigen::MatrixXd > > & potentialParents,
                          const Eigen::VectorXd & markerImbalances, const double & convexMarkerImbalanceInterpolation, const Eigen::VectorXd & tolerance, const double & theta, const Eigen::VectorXd & alleleFrequencies,
                          const std::size_t & levelsOfStutterRecursion);

        Eigen::VectorXd GenerateUnknownGenotype(const std::size_t & seed);

};

#endif

