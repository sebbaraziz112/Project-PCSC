#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <gtest/gtest.h>

#include "../includes/utils.hpp"
#include "../includes/methods.hpp"



TEST(MATRIX_VECTOR_TESTS, FACTORS){
    
    const int N = 183839;
    std::shared_ptr<AbstractMethod> method = std::make_shared<IdentityMethod>();
    std::vector<int> factors = method->getFactors(N);
    EXPECT_EQ(factors[0] * factors[1], N);

}

TEST(MATRIX_VECTOR_TESTS, IDENTITY){

    const int N = 183839;
    std::shared_ptr<AbstractMethod> method = std::make_shared<IdentityMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> vector_matrices;
    
    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> mat1(1, N);
    randomInit<uint16_t>(mat1);
    vector_matrices.push_back(mat1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> mat2(1, N);
    randomInit<uint16_t>(mat2);
    vector_matrices.push_back(mat2);

    std::vector<int> factors = method->getFactors(N);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> matrix_matrices = method->vectorToMatrix<uint16_t>(vector_matrices, factors);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> result_vector_matrices = method->matrixToVector<uint16_t>(matrix_matrices);

    EXPECT_EQ(vector_matrices[0], result_vector_matrices[0]);

}

