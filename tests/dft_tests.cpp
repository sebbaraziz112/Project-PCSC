#include "../includes/methods.hpp"
#include "../includes/utils.hpp"
#include "../includes/DataReader.hpp"
#include <gtest/gtest.h>


// Radix-2 Method Tests

TEST(FFT1D_TEST, FORWARD_PASS_1D){

    std::shared_ptr<FFT1DMethod> myFFT = std::make_shared<FFT1DMethod>(false);
    std::shared_ptr<NaiveDFTMethod> myDFT = std::make_shared<NaiveDFTMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 16);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 16);
    randomInit(test2);
    tests.push_back(test2);

    

    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> tests_truth = myDFT->computeMethod<uint16_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultFFT = myFFT->computeMethod<uint16_t, std::complex<double>>(tests);


    for (int j = 0; j < 16; j++) {
        EXPECT_NEAR(resultFFT[0](0, j).real(), tests_truth[0](0, j).real(), 1e-3);
        EXPECT_NEAR(resultFFT[0](0, j).imag(), tests_truth[0](0, j).imag(), 1e-3);
    }

    for (int j = 0; j < 11; j++) {
        EXPECT_NEAR(resultFFT[1](0, j).real(), tests_truth[1](0, j).real(), 1e-3);
        EXPECT_NEAR(resultFFT[1](0, j).imag(), tests_truth[1](0, j).imag(), 1e-3);
    }

}

TEST(FFT1D_TEST, BACKWARD_PASS_1D){

    std::shared_ptr<Inv1DFFTMethod> myInvFFT = std::make_shared<Inv1DFFTMethod>(false);
    std::shared_ptr<NaiveDFTMethod> myDFT = std::make_shared<NaiveDFTMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 16);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 16);
    randomInit(test2);
    tests.push_back(test2);

    

    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> tests_inverse_dft = myDFT->computeMethod<uint16_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvFFT = myInvFFT->computeMethod<std::complex<double>, uint16_t>(tests_inverse_dft);


    EXPECT_EQ(resultInvFFT[0], tests[0]);
    EXPECT_EQ(resultInvFFT[1], tests[1]);

}

TEST(FFT1D_TEST, IDENTITY_1D){
    std::shared_ptr<FFT1DMethod> myFFT = std::make_shared<FFT1DMethod>(false);
    std::shared_ptr<Inv1DFFTMethod> myInvFFT = std::make_shared<Inv1DFFTMethod>(false);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 16);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 16);
    randomInit(test2);
    tests.push_back(test2);


    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultFFT = myFFT->computeMethod<uint16_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> invFFT = myInvFFT->computeMethod<std::complex<double>, uint16_t>(resultFFT);


    EXPECT_EQ(invFFT[0], tests[0]);
    EXPECT_EQ(invFFT[1], tests[1]);

}

////////////////////////////////////////

// Bluestein1D Method Test

TEST(BLUESTEIN1D_TEST, FORWARD_PASS_1D){

    std::shared_ptr<BlueStein1DMethod> myFFT = std::make_shared<BlueStein1DMethod>();
    std::shared_ptr<NaiveDFTMethod> myDFT = std::make_shared<NaiveDFTMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 17);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 17);
    randomInit(test2);
    tests.push_back(test2);

    

    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> tests_truth = myDFT->computeMethod<uint16_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultFFT = myFFT->computeMethod<uint16_t, std::complex<double>>(tests);


    for (int j = 0; j < 16; j++) {
        EXPECT_NEAR(resultFFT[0](0, j).real(), tests_truth[0](0, j).real(), 1e-3);
        EXPECT_NEAR(resultFFT[0](0, j).imag(), tests_truth[0](0, j).imag(), 1e-3);
    }

    for (int j = 0; j < 11; j++) {
        EXPECT_NEAR(resultFFT[1](0, j).real(), tests_truth[1](0, j).real(), 1e-3);
        EXPECT_NEAR(resultFFT[1](0, j).imag(), tests_truth[1](0, j).imag(), 1e-3);
    }

}

TEST(BLUESTEIN1D_TEST, BACKWARD_PASS_1D){

    std::shared_ptr<InvBlueStein1DMethod> myInvFFT = std::make_shared<InvBlueStein1DMethod>();
    std::shared_ptr<NaiveDFTMethod> myDFT = std::make_shared<NaiveDFTMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 17);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 17);
    randomInit(test2);
    tests.push_back(test2);

    

    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> tests_inverse_dft = myDFT->computeMethod<uint16_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvFFT = myInvFFT->computeMethod<std::complex<double>, uint16_t>(tests_inverse_dft);


    EXPECT_EQ(resultInvFFT[0], tests[0]);
    EXPECT_EQ(resultInvFFT[1], tests[1]);

}

TEST(BLUESTEIN1D_TEST, IDENTITY_1D){
    std::shared_ptr<BlueStein1DMethod> myFFT = std::make_shared<BlueStein1DMethod>();
    std::shared_ptr<InvBlueStein1DMethod> myInvFFT = std::make_shared<InvBlueStein1DMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 17);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 17);
    randomInit(test2);
    tests.push_back(test2);


    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultFFT = myFFT->computeMethod<uint16_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> invFFT = myInvFFT->computeMethod<std::complex<double>, uint16_t>(resultFFT);


    EXPECT_EQ(invFFT[0], tests[0]);
    EXPECT_EQ(invFFT[1], tests[1]);

}

/////////////////////////////////////////

// Bluestein Method Test

TEST(BLUESTEIN_TEST, FORWARD_PASS_1D){

    std::shared_ptr<BlueSteinMethod> myBlueStein = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<NaiveDFTMethod> myDFT = std::make_shared<NaiveDFTMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 11);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 11);
    randomInit(test2);
    tests.push_back(test2);

    

    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> tests_truth = myDFT->computeMethod<uint16_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBlueStein = myBlueStein->computeMethod<uint16_t, std::complex<double>>(tests);


    for (int j = 0; j < 11; j++) {
        EXPECT_NEAR(resultBlueStein[0](0, j).real(), tests_truth[0](0, j).real(), 1e-3);
        EXPECT_NEAR(resultBlueStein[0](0, j).imag(), tests_truth[0](0, j).imag(), 1e-3);
    }

    for (int j = 0; j < 11; j++) {
        EXPECT_NEAR(resultBlueStein[1](0, j).real(), tests_truth[1](0, j).real(), 1e-3);
        EXPECT_NEAR(resultBlueStein[1](0, j).imag(), tests_truth[1](0, j).imag(), 1e-3);
    }

}

TEST(BLUESTEIN_TEST, FORWARD_PASS_2D){

    std::shared_ptr<BlueSteinMethod> myBlueStein = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<NaiveDFTMethod> myDFT = std::make_shared<NaiveDFTMethod>();

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(11, 11);
    randomInit<uint8_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(11, 11);
    randomInit<uint8_t>(test2);
    tests.push_back(test2);

    

    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> tests_truth = myDFT->computeMethod<uint8_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBlueStein = myBlueStein->computeMethod<uint8_t, std::complex<double>>(tests);


    for (int j = 0; j < 11; j++) {
        for (int i = 0; i < 11; i ++){
            EXPECT_NEAR(resultBlueStein[0](i, j).real(), tests_truth[0](i, j).real(), 1e-3);
            EXPECT_NEAR(resultBlueStein[0](i, j).imag(), tests_truth[0](i, j).imag(), 1e-3);
        }
    }

    for (int j = 0; j < 11; j++) {
        for (int i = 0; i < 11; i ++){
            EXPECT_NEAR(resultBlueStein[1](i, j).real(), tests_truth[1](i, j).real(), 1e-3);
            EXPECT_NEAR(resultBlueStein[1](i, j).imag(), tests_truth[1](i, j).imag(), 1e-3);
        }
    }

}

TEST(BLUESTEIN_TEST, BACKWARD_PASS_1D){

    std::shared_ptr<InvBlueSteinMethod> myInvFFT = std::make_shared<InvBlueSteinMethod>();
    std::shared_ptr<NaiveDFTMethod> myDFT = std::make_shared<NaiveDFTMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(1, 11);
    randomInit<uint16_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(1, 11);
    randomInit(test2);
    tests.push_back(test2);

    

    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> tests_inverse_dft = myDFT->computeMethod<uint16_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvFFT = myInvFFT->computeMethod<std::complex<double>, uint16_t>(tests_inverse_dft);


    EXPECT_EQ(resultInvFFT[0], tests[0]);
    EXPECT_EQ(resultInvFFT[1], tests[1]);

}

TEST(BLUESTEIN_TEST, BACKWARD_PASS_2D){

    std::shared_ptr<InvBlueSteinMethod> myInvFFT = std::make_shared<InvBlueSteinMethod>();
    std::shared_ptr<NaiveDFTMethod> myDFT = std::make_shared<NaiveDFTMethod>();

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> tests;

    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> test1;
    test1.resize(11, 11);
    randomInit<uint8_t>(test1);
    tests.push_back(test1);

    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> test2;
    test2.resize(11, 11);
    randomInit(test2);
    tests.push_back(test2);

    

    
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> tests_inverse_dft = myDFT->computeMethod<uint8_t, std::complex<double>>(tests);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvFFT = myInvFFT->computeMethod<std::complex<double>, uint8_t>(tests_inverse_dft);


    EXPECT_EQ(resultInvFFT[0], tests[0]);
    EXPECT_EQ(resultInvFFT[1], tests[1]);

}

TEST(BLUESTEIN_TEST, IDENTITY_1D){

    std::shared_ptr<BlueSteinMethod> myBlueStein = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<InvBlueSteinMethod> myInvBlueStein = std::make_shared<InvBlueSteinMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests_1d;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1_1d;
    test1_1d.resize(1, 10);
    randomInit<uint16_t>(test1_1d);
    tests_1d.push_back(test1_1d);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2_1d;
    test2_1d.resize(1, 10);
    randomInit<uint16_t>(test2_1d);
    tests_1d.push_back(test2_1d);

    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBluestein_1d = myBlueStein->computeMethod<uint16_t, std::complex<double>>(tests_1d);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_1d = myInvBlueStein->computeMethod<std::complex<double>, uint16_t>(resultBluestein_1d);


    EXPECT_EQ(resultInvBluestein_1d[0], tests_1d[0]);
    EXPECT_EQ(resultInvBluestein_1d[1], tests_1d[1]);
}

TEST(BLUESTEIN_TEST, IDENTITY_2D){
    std::shared_ptr<BlueSteinMethod> myBlueStein = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<InvBlueSteinMethod> myInvBlueStein = std::make_shared<InvBlueSteinMethod>();

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> tests_2d;

    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> test1_2d;
    test1_2d.resize(10, 10);
    randomInit<uint8_t>(test1_2d);
    tests_2d.push_back(test1_2d);

    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> test2_2d;
    test2_2d.resize(10, 10);
    randomInit<uint8_t>(test2_2d);
    tests_2d.push_back(test2_2d);

    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBluestein_2d = myBlueStein->computeMethod<uint8_t, std::complex<double>>(tests_2d);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_2d = myInvBlueStein->computeMethod<std::complex<double>, uint8_t>(resultBluestein_2d);

    EXPECT_EQ(resultInvBluestein_2d[0], tests_2d[0]);
    EXPECT_EQ(resultInvBluestein_2d[1], tests_2d[1]);
}

TEST(BLUESTEIN_TEST, IDENTITY_1D_BIG){

    std::shared_ptr<BlueSteinMethod> myBlueStein = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<InvBlueSteinMethod> myInvBlueStein = std::make_shared<InvBlueSteinMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests_1d;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1_1d;
    test1_1d.resize(1, 1000);
    randomInit<uint16_t>(test1_1d);
    tests_1d.push_back(test1_1d);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2_1d;
    test2_1d.resize(1, 1000);
    randomInit<uint16_t>(test2_1d);
    tests_1d.push_back(test2_1d);

    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBluestein_1d = myBlueStein->computeMethod<uint16_t, std::complex<double>>(tests_1d);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_1d = myInvBlueStein->computeMethod<std::complex<double>, uint16_t>(resultBluestein_1d);


    EXPECT_EQ(resultInvBluestein_1d[0], tests_1d[0]);
    EXPECT_EQ(resultInvBluestein_1d[1], tests_1d[1]);
}

TEST(BLUESTEIN_TEST, IDENTITY_2D_BIG){
    std::shared_ptr<BlueSteinMethod> myBlueStein = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<InvBlueSteinMethod> myInvBlueStein = std::make_shared<InvBlueSteinMethod>();

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> tests_2d;

    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> test1_2d;
    test1_2d.resize(100, 100);
    randomInit<uint8_t>(test1_2d);
    tests_2d.push_back(test1_2d);

    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> test2_2d;
    test2_2d.resize(100, 100);
    randomInit<uint8_t>(test2_2d);
    tests_2d.push_back(test2_2d);

    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBluestein_2d = myBlueStein->computeMethod<uint8_t, std::complex<double>>(tests_2d);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_2d = myInvBlueStein->computeMethod<std::complex<double>, uint8_t>(resultBluestein_2d);

    EXPECT_EQ(resultInvBluestein_2d[0], tests_2d[0]);
    EXPECT_EQ(resultInvBluestein_2d[1], tests_2d[1]);
}

TEST(LINE_ONLY_BLUESTEIN_TEST, IDENTITY_1D){

    std::shared_ptr<LineOnlyBlueSteinMethod> myBlueStein = std::make_shared<LineOnlyBlueSteinMethod>();
    std::shared_ptr<LineOnlyInvBlueSteinMethod> myInvBlueStein = std::make_shared<LineOnlyInvBlueSteinMethod>();

    const int N = 100;
    std::vector<int> factors = myBlueStein->getFactors(N);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests_2d;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1_2d;
    test1_2d.resize(1, N);
    randomInit<uint16_t>(test1_2d);
    tests_2d.push_back(test1_2d);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2_2d;
    test2_2d.resize(1, N);
    randomInit<uint16_t>(test2_2d);
    tests_2d.push_back(test2_2d);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests_2d_matrices = myBlueStein->vectorToMatrix<uint16_t>(tests_2d, factors);


    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBluestein_2d = myBlueStein->computeMethod<uint16_t, std::complex<double>>(tests_2d_matrices);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_2d = myInvBlueStein->computeMethod<std::complex<double>, uint16_t>(resultBluestein_2d);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_2d_vectors = myBlueStein->matrixToVector<uint16_t>(resultInvBluestein_2d);

    EXPECT_EQ(resultInvBluestein_2d_vectors[0], tests_2d[0]);
    EXPECT_EQ(resultInvBluestein_2d_vectors[1], tests_2d[1]);

}

TEST(LINE_ONLY_BLUESTEIN_TEST, IDENTITY_1D_BIG){

    std::shared_ptr<LineOnlyBlueSteinMethod> myBlueStein = std::make_shared<LineOnlyBlueSteinMethod>();
    std::shared_ptr<LineOnlyInvBlueSteinMethod> myInvBlueStein = std::make_shared<LineOnlyInvBlueSteinMethod>();

    const int N = 183839;
    std::vector<int> factors = myBlueStein->getFactors(N);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests_2d;

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test1_2d;
    test1_2d.resize(1, N);
    randomInit16<uint16_t>(test1_2d);
    tests_2d.push_back(test1_2d);

    Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> test2_2d;
    test2_2d.resize(1, N);
    randomInit16<uint16_t>(test2_2d);
    tests_2d.push_back(test2_2d);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests_2d_matrices = myBlueStein->vectorToMatrix<uint16_t>(tests_2d, factors);


    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBluestein_2d = myBlueStein->computeMethod<uint16_t, std::complex<double>>(tests_2d_matrices);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_2d = myInvBlueStein->computeMethod<std::complex<double>, uint16_t>(resultBluestein_2d);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_2d_vectors = myBlueStein->matrixToVector<uint16_t>(resultInvBluestein_2d);

    EXPECT_EQ(resultInvBluestein_2d_vectors[0], tests_2d[0]);
    EXPECT_EQ(resultInvBluestein_2d_vectors[1], tests_2d[1]);

}

TEST(LINE_ONLY_BLUESTEIN_TEST, IDENTITY_1D_SOUND){

    std::shared_ptr<LineOnlyBlueSteinMethod> myBlueStein = std::make_shared<LineOnlyBlueSteinMethod>();
    std::shared_ptr<LineOnlyInvBlueSteinMethod> myInvBlueStein = std::make_shared<LineOnlyInvBlueSteinMethod>();
    std::shared_ptr<SoundReader> mySoundReader = std::make_shared<SoundReader>(true, true);

    

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests_2d = mySoundReader->getData("../ressources/TrackFolder/noisy1.wav");

    const int N = tests_2d[0].cols();
    std::vector<int> factors = myBlueStein->getFactors(N);

    
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> tests_2d_matrices = myBlueStein->vectorToMatrix<uint16_t>(tests_2d, factors);


    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> resultBluestein_2d = myBlueStein->computeMethod<uint16_t, std::complex<double>>(tests_2d_matrices);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_2d = myInvBlueStein->computeMethod<std::complex<double>, uint16_t>(resultBluestein_2d);

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> resultInvBluestein_2d_vectors = myBlueStein->matrixToVector<uint16_t>(resultInvBluestein_2d);

    EXPECT_EQ(resultInvBluestein_2d_vectors[0], tests_2d[0]);

}

/////////////////////////////////////////