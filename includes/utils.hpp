/**
 * @file MatrixUtils.hpp
 * @brief Utility functions for matrix operations, including random initialization,
 *        type casting, file comparison, and display preparation.
 *
 * This file provides helper functions and templates to handle Eigen matrices,
 * initialize them with random values, cast between types, and compare files.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <cmath>

/**
 * @brief Compare the contents of two files line by line.
 * @param file1_name Path to the first file.
 * @param file2_name Path to the second file.
 * @return true if the files are identical, false otherwise.
 */
bool compareFiles(std::string file1_name, std::string file2_name);

/**
 * @brief Display or prepare matrices for visualization or processing.
 * @param upperBound Maximum value for display scaling or normalization.
 * @param matrices Vector of Eigen matrices to prepare.
 *
 * This function can be used to normalize, scale, or display matrices
 * for debugging or visualization purposes.
 */
void displayPreparation(double upperBound, std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>& matrices);

/**
 * @brief Initialize a matrix with random integer values in the range [0, 255].
 * @tparam T Type of the matrix elements (e.g., int, double, float).
 * @param matrix Eigen matrix to initialize.
 *
 * Uses a uniform random distribution and fills each element with
 * a random value cast to type T.
 */
template <typename T>
void randomInit(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix){
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    for (int i = 0; i < nRows; i++){
        for (int j = 0; j < nCols; j++){
            matrix(i, j) = static_cast<T>(dist(gen));
        }
    }
}

/**
 * @brief Cast a vector of matrices from one type to another.
 * @tparam I Input matrix element type.
 * @tparam O Output matrix element type.
 * @param matrices Vector of input matrices of type I.
 * @return Vector of matrices of type O with the same dimensions.
 *
 * If the input type I and output type O are identical, returns a copy of the input vector.
 * Otherwise, each element is cast to type O, and the absolute value is rounded
 * before casting.
 */
template <typename I, typename O>
std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> CastMatrix(
    std::vector<Eigen::Matrix<I, Eigen::Dynamic, Eigen::Dynamic>>& matrices) 
{
    if constexpr (std::is_same_v<O, I>) {
        std::vector<Eigen::Matrix<I, Eigen::Dynamic, Eigen::Dynamic>> OVector;
        for (auto& mat : matrices) {
            OVector.push_back(mat);  
        }
        return OVector;
    }

    std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> OVector;
    const int nRows = matrices[0].rows();
    const int nCols = matrices[0].cols();
    
    for (int i = 0; i < matrices.size(); i++){
        Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic> mat(nRows, nCols);
        OVector.push_back(mat);
    }

    for (int chann = 0; chann < matrices.size(); chann++){
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                OVector[chann](i, j) = static_cast<O>(std::round(std::abs(matrices[chann](i, j))));
            }
        }
    }

    return OVector;
}
