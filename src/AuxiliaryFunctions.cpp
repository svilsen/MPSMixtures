#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/math/special_functions/binomial.hpp>

//[[Rcpp::export(.partialSumEigen)]]
Eigen::VectorXd partialSumEigen(const Eigen::VectorXd & x)
{
    const std::size_t N = x.size();
    Eigen::VectorXd partialSum(N + 1);

    partialSum[0] = 0;
    for (std::size_t i = 1; i < N + 1; i++)
    {
        partialSum[i] = x[i - 1] + partialSum[i - 1];
    }

    return partialSum;
}

Eigen::MatrixXd bindColumns(const Eigen::MatrixXd & A, const Eigen::MatrixXd & B)
{
    Eigen::MatrixXd C(A.rows(), A.cols() + B.cols());
    C << A, B;
    return C;
}


Eigen::VectorXd STDEigen(const std::vector<double> &x)
{
    Eigen::VectorXd xEigen(x.size());
    for (std::size_t i = 0; i < x.size(); i++)
        xEigen[i] = x[i];

    return xEigen;
}

std::vector<double> EigenSTD(const Eigen::VectorXd &x)
{
    std::vector<double> xNEWMAT(x.size());
    for (std::size_t i = 0; i < x.size(); i++)
        xNEWMAT[i] = x[i];

    return xNEWMAT;
}

Eigen::Vector2i nonZeroElementsOfMarker(const Eigen::VectorXd & decodedProfile_m)
{
    const std::size_t N = decodedProfile_m.size();

    Eigen::Vector2i res;
    std::size_t i = 0;
    std::size_t j = 0;
    while(j < 2)
    {
        if (decodedProfile_m[i] != 0)
        {
            int k = static_cast<int>(decodedProfile_m[i]);
            while(k > 0)
            {
                res[j] = i;
                k -= 1;
                j += 1;
            }
        }
        i += 1;
    }

    return res;
}

std::vector<int> sortedIndex(const Eigen::VectorXd & x)
{
    std::vector<int> x_sorted(x.size());
    std::iota(x_sorted.begin(), x_sorted.end(), 0);
    auto comparator = [&x](int i, int j){ return x[i] > x[j]; };

    std::sort(x_sorted.begin(), x_sorted.end(), comparator);

    return x_sorted;
}

Eigen::VectorXd noiseQuantiles(const Eigen::VectorXd & Coverage,
                               const Eigen::VectorXd & PartialSumAlleles,
                               const std::vector<Eigen::VectorXd> & NoiseIndex,
                               const int & NoiseSize, const double & q)
{
    int k = 0;
    Eigen::VectorXd noiseSorted(NoiseSize);
    for (std::size_t m = 0; m < NoiseIndex.size(); m++)
    {
        const Eigen::VectorXd & NoiseIndex_m = NoiseIndex[m];
        for (std::size_t a = 0; a < NoiseIndex_m.size(); a++)
        {
            std::size_t n = PartialSumAlleles[m] + NoiseIndex_m[a];

            int h = k - 1;
            while ((h >= 0) & (noiseSorted[h] > Coverage[n])) {
                noiseSorted[h + 1] = noiseSorted[h];
                h--;
            }

            noiseSorted[h + 1] = Coverage[n];
            k++;
        }
    }

    const double a = (1.0 - q) / 2.0;
    int li = std::floor(a * NoiseSize);
    int ui = std::ceil((1 - a) * NoiseSize);
    if (ui > (NoiseSize - 1)) {
        ui = NoiseSize - 1;
    }

    Eigen::VectorXd quantiles(2);
    quantiles[0] = noiseSorted[li];
    quantiles[1] = noiseSorted[ui];
    return quantiles;
}

//[[Rcpp::export()]]
Eigen::MatrixXd generatePossibleGenotypes(const std::size_t & N)
{
    if (N == 1)
    {
        Eigen::MatrixXd possibleSingleprofilesMatrix = 2*Eigen::MatrixXd::Ones(N, N);
        return possibleSingleprofilesMatrix;
    }

    const std::size_t possibleSingles = N;
    const std::size_t possiblePairs = boost::math::binomial_coefficient<double>(N, 2);

    const std::size_t totalPossibleOutcomes = possibleSingles + possiblePairs;
    const std::size_t totalRuns = std::floor(possiblePairs / N);

    Eigen::MatrixXd possibleSingleprofilesMatrix = Eigen::MatrixXd::Zero(N, totalPossibleOutcomes);
    for (std::size_t i = 0; i < N; i++)
    {
        // Single Values
        Eigen::VectorXd singleVector = Eigen::VectorXd::Zero(N);
        singleVector(i) = 2;
        possibleSingleprofilesMatrix.col(i) = singleVector;

        // Paired Values
        for (std::size_t j = 1; j <= totalRuns; j++)
        {
            Eigen::VectorXd pairsVector = Eigen::VectorXd::Zero(N);
            pairsVector(i) = 1;
            pairsVector((i + j) % N) = 1;
            possibleSingleprofilesMatrix.col(j*N + i) = pairsVector;
        }
    }

    // Remaining Paired Values
    for (std::size_t k = (totalRuns + 1)*N; k < totalPossibleOutcomes; k++)
    {
        Eigen::VectorXd pairsVector = Eigen::VectorXd::Zero(N);
        pairsVector(k % N) = 1;
        pairsVector(((k % N) + (totalRuns + 1)) % N) = 1;
        possibleSingleprofilesMatrix.col(k) = pairsVector;
    }

    return possibleSingleprofilesMatrix;
}
