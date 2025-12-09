#include "../includes/methods.hpp"
#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/utils.hpp"
#include <gtest/gtest.h>
#include <fstream>


TEST(BLUESTEIN_TEST_IMAGE, IDENTITY){
    std::shared_ptr<BlueSteinMethod> myBlueStein = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<InvBlueSteinMethod> myInvBlueStein = std::make_shared<InvBlueSteinMethod>();

    std::shared_ptr<ImageReader> myImageReader = std::make_shared<ImageReader>(true, true);
    std::shared_ptr<ImageWriter> myImageWriter = std::make_shared<ImageWriter>();

    std::string original_file_image = "../ressources/ImageFolder/Mauritius_beach.png";

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> tests_2d = myImageReader->getData(original_file_image);

    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBluestein_2d = myBlueStein->computeMethod<uint8_t, std::complex<double>>(tests_2d);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_2d = myInvBlueStein->computeMethod<std::complex<double>, uint8_t>(resultBluestein_2d);

    EXPECT_EQ(resultInvBluestein_2d[0], tests_2d[0]);
    EXPECT_EQ(resultInvBluestein_2d[1], tests_2d[1]);
    EXPECT_EQ(resultInvBluestein_2d[2], tests_2d[2]);

    
}
