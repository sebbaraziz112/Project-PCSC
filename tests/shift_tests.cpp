#include "../includes/methods.hpp"
#include "../includes/utils.hpp"
#include <gtest/gtest.h>



TEST(SHIFT_MATRIX, IDENTITY_EVEN_TEST_2D){

    std::shared_ptr<BandFiltering> myBandFilter = std::make_shared<BandFiltering>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(16, 16);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(16, 16);
    randomInit(test2);
    tests.push_back(test2);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> shifted = myBandFilter->shiftCenter<uint16_t>(tests);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> invShifted = myBandFilter->shiftInverse<uint16_t>(shifted);

    EXPECT_EQ(invShifted[0], tests[0]);
    EXPECT_EQ(invShifted[1], tests[1]);

}

TEST(SHIFT_MATRIX, IDENTITY_ODD_TEST_2D){

    std::shared_ptr<BandFiltering> myBandFilter = std::make_shared<BandFiltering>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(17, 17);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 17);
    randomInit(test2);
    tests.push_back(test2);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> shifted = myBandFilter->shiftCenter<uint16_t>(tests);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> invShifted = myBandFilter->shiftInverse<uint16_t>(shifted);

    EXPECT_EQ(invShifted[0], tests[0]);
    EXPECT_EQ(invShifted[1], tests[1]);
    
}

TEST(SHIFT_MATRIX, IDENTITY_EVEN_TEST_1D){

    std::shared_ptr<BandFiltering> myBandFilter = std::make_shared<BandFiltering>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 16);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 16);
    randomInit(test2);
    tests.push_back(test2);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> shifted = myBandFilter->shiftCenter<uint16_t>(tests);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> invShifted = myBandFilter->shiftInverse<uint16_t>(shifted);

    EXPECT_EQ(invShifted[0], tests[0]);
    EXPECT_EQ(invShifted[1], tests[1]);

}

TEST(SHIFT_MATRIX, IDENTITY_ODD_TEST_1D){

    std::shared_ptr<BandFiltering> myBandFilter = std::make_shared<BandFiltering>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 17);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 17);
    randomInit(test2);
    tests.push_back(test2);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> shifted = myBandFilter->shiftCenter<uint16_t>(tests);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> invShifted = myBandFilter->shiftInverse<uint16_t>(shifted);

    EXPECT_EQ(invShifted[0], tests[0]);
    EXPECT_EQ(invShifted[1], tests[1]);
    
}





