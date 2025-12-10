#include <iostream>
#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/methods.hpp"
#include "../includes/utils.hpp"
#include "../includes/DataTypes.hpp"

int main(int argc, char* argv[]){
    if (argc < 4) {
        std::cout << "Usage: imsonpro bandfilter-sound <file> <prct bottom> <prct top>\n";
        return 0;
    }

    std::string input_filename = argv[1];
    double prct_bottom = std::stod(argv[2]);
    double prct_top = std::stod(argv[3]);

    std::shared_ptr<SoundReader> mySoundReader = std::make_shared<SoundReader>(true, true);
    std::shared_ptr<SoundWriter> mySoundWriter = std::make_shared<SoundWriter>();
    std::shared_ptr<BandFiltering> myBandFilter = std::make_shared<BandFiltering>();
    

    myBandFilter->setLowerFractionX(0.0);
    myBandFilter->setUpperFractionX(1.0);
    myBandFilter->setLowerFractionY(prct_bottom / 100.0);
    myBandFilter->setUpperFractionY(prct_top / 100.0);
    myBandFilter->setAmplificationFactor(1.0);

    

    std::string output_filename = mySoundReader->getPreExtension(input_filename) + "_bandfiltered.wav";

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundData = mySoundReader->getData(input_filename);

    if (soundData[0].cols() < 10000){

        std::shared_ptr<BlueSteinMethod> myBS = std::make_shared<BlueSteinMethod>();
        std::shared_ptr<InvBlueSteinMethod> myIBS = std::make_shared<InvBlueSteinMethod>();
        std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> soundDataDFT = myBS->computeMethod<uint16_t, std::complex<double>>(soundData);
        std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> shiftedDFT = myBandFilter->shiftCenter<std::complex<double>>(soundDataDFT);
        std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> bandFilteredDFT = myBandFilter->computeMethod<std::complex<double>, std::complex<double>>(shiftedDFT);
        std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> unshiftedDFT = myBandFilter->shiftInverse<std::complex<double>>(bandFilteredDFT);
        std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundDataIDFT = myIBS->computeMethod<std::complex<double>, uint16_t>(unshiftedDFT);
        
        WAVHeader myHeaderSound = mySoundReader->getFinalHeader(input_filename);
        mySoundWriter->write(soundDataIDFT, output_filename, myHeaderSound);

    } else {

        const int N = soundData[0].cols();
        std::shared_ptr<ImageWriter> myImageWriter = std::make_shared<ImageWriter>();

        std::string dft_image_filename = mySoundReader->getPreExtension(input_filename) + "_dft_image.bmp";
        std::string dft_banded_image_filename = mySoundReader->getPreExtension(input_filename) + "_dft_banded_image.bmp";

        std::shared_ptr<LineOnlyBlueSteinMethod> myBS = std::make_shared<LineOnlyBlueSteinMethod>();
        std::shared_ptr<LineOnlyInvBlueSteinMethod> myIBS = std::make_shared<LineOnlyInvBlueSteinMethod>();

        std::vector<int> factors = myBS->getFactors(N);

        std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> matrixData = myBS->vectorToMatrix<uint16_t>(soundData, factors);
        std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> soundDataDFT = myBS->computeMethod<uint16_t, std::complex<double>>(matrixData);
        std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> shiftedDFT = myBandFilter->shiftCenter<std::complex<double>>(soundDataDFT);
        std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> bandFilteredDFT = myBandFilter->computeMethod<std::complex<double>, std::complex<double>>(shiftedDFT);
        std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> unshiftedDFT = myBandFilter->shiftInverse<std::complex<double>>(bandFilteredDFT);
        std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundDataIDFT = myIBS->computeMethod<std::complex<double>, uint16_t>(unshiftedDFT);
        std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundDataIDFTVector = myBS->matrixToVector<uint16_t>(soundDataIDFT);

        WAVHeader myHeaderSound = mySoundReader->getFinalHeader(input_filename);
        mySoundWriter->write(soundDataIDFTVector, output_filename, myHeaderSound);

        std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> dftImageData = CastMatrix<std::complex<double>, double>(shiftedDFT);
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> dftBandedImageData = CastMatrix<std::complex<double>, double>(bandFilteredDFT);
        displayPreparation(255.0, dftImageData);
        displayPreparation(255.0, dftBandedImageData);
        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> dftFinalImageData = CastMatrix<double, uint8_t>(dftImageData);
        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> dftBandedFinalImageData = CastMatrix<double, uint8_t>(dftBandedImageData);

        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> tripleddftFinalImageData(3, dftFinalImageData[0]);
        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> tripleddftBandedFinalImageData(3, dftBandedFinalImageData[0]);


        BMPHeader myHeaderImage;
        initBMPHeader(myHeaderImage, dftFinalImageData[0].cols(), dftFinalImageData[0].rows(), 24);
        
        myImageWriter->write(tripleddftFinalImageData, dft_image_filename, myHeaderImage);
        myImageWriter->write(tripleddftBandedFinalImageData, dft_banded_image_filename, myHeaderImage);

    }  
    
}