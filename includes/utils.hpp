#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <random>


bool compareFiles(std::string file1_name, std::string file2_name);

void displayPreparation(double upperBound, std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>& matrices);

template <typename T>
void randomInit(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix){
    const int nRows = matrix.rows();
    const int nCols = matrix.cols();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);


    for (int i = 0; i < nRows; i ++){
        for (int j = 0; j < nCols; j ++){
            matrix(i, j) = static_cast<T>(dist(gen));
        }
    }
}

template <typename I, typename O>
std::vector<Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic>> CastMatrix(std::vector<Eigen::Matrix<I, Eigen::Dynamic, Eigen::Dynamic>>& matrices) {
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
    for (int i = 0; i < matrices.size(); i ++){
        Eigen::Matrix<O, Eigen::Dynamic, Eigen::Dynamic> mat(nRows, nCols);
        OVector.push_back(mat);
    }
    for (int chann = 0; chann < matrices.size(); chann ++){
        for (int i = 0; i < nRows; i ++){
            for (int j = 0; j < nCols; j ++){
                OVector[chann](i, j) = static_cast<O>(std::round(std::abs(matrices[chann](i, j))));
            }
        }
    }
    return OVector;
}


template <typename T>
std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> shiftCenter(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) {
    const int nChann = matrices.size();
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> centeredMatrices;
    centeredMatrices.reserve(nChann);

    for (int chann = 0; chann < nChann; chann++) {
        const auto& mat = matrices[chann];
        const int nRows = mat.rows();
        const int nCols = mat.cols();

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> newMatrix(nRows, nCols);

        int midRow = nRows / 2;
        int midCol = nCols / 2;
        newMatrix.topLeftCorner(midRow, midCol) = mat.bottomRightCorner(nRows - midRow, nCols - midCol);
        newMatrix.topRightCorner(midRow, nCols - midCol) = mat.bottomLeftCorner(nRows - midRow, midCol);
        newMatrix.bottomLeftCorner(nRows - midRow, midCol) = mat.topRightCorner(midRow, nCols - midCol);
        newMatrix.bottomRightCorner(nRows - midRow, nCols - midCol) = mat.topLeftCorner(midRow, midCol);

        centeredMatrices.push_back(newMatrix);
    }

    return centeredMatrices;
}


template <typename T>
std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> shiftInverse(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrices) 
{
    const int nChann = matrices.size();
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> invMatrices;
    invMatrices.reserve(nChann);

    for (int chann = 0; chann < nChann; chann++) {
        const auto& mat = matrices[chann];
        const int nRows = mat.rows();
        const int nCols = mat.cols();

        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> newMatrix(nRows, nCols);

        int midRow = nRows / 2;
        int midCol = nCols / 2;

        newMatrix.topLeftCorner(midRow, midCol) = mat.bottomRightCorner(nRows - midRow, nCols - midCol);
        newMatrix.topRightCorner(midRow, nCols - midCol) = mat.bottomLeftCorner(nRows - midRow, midCol);
        newMatrix.bottomLeftCorner(nRows - midRow, midCol) = mat.topRightCorner(midRow, nCols - midCol);
        newMatrix.bottomRightCorner(nRows - midRow, nCols - midCol) = mat.topLeftCorner(midRow, midCol);

        invMatrices.push_back(newMatrix);
    }

    return invMatrices;
}

