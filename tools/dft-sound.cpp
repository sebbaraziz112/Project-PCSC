#include <iostream>
#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/methods.hpp"
#include "../includes/utils.hpp"
#include "../includes/DataTypes.hpp"

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Usage: imsonpro dft-sound <file> \n";
        return 0;
    }

    std::string input_filename = argv[1];

    std::shared_ptr<SoundReader> mySoundReader = std::make_shared<SoundReader>(true, true);
    std::shared_ptr<SoundWriter> mySoundWriter = std::make_shared<SoundWriter>();

    

    std::string output_filename = mySoundReader->getPreExtension(input_filename) + "_dft_image.bmp";

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundData = mySoundReader->getData(input_filename);
    std::shared_ptr<ImageWriter> myImageWriter = std::make_shared<ImageWriter>();
    std::shared_ptr<BandFiltering> myBandFilter = std::make_shared<BandFiltering>();


    if (soundData[0].cols() < 10000){

        std::shared_ptr<BlueSteinMethod> myBS = std::make_shared<BlueSteinMethod>();
        std::shared_ptr<InvBlueSteinMethod> myIBS = std::make_shared<InvBlueSteinMethod>();
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> soundDataDFT = myBS->computeMethod<uint16_t, double>(soundData);
        displayPreparation(255.0, soundDataDFT);
        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> soundDataDFTUint8 = CastMatrix<double, uint8_t>(soundDataDFT);
        
        BMPHeader myHeaderImage;
        initBMPHeader(myHeaderImage, soundDataDFT[0].cols(), soundDataDFT[0].rows(), 24);

        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> tripleddftFinalImageData(3, soundDataDFTUint8[0]);
        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> finalDFT = myBandFilter->shiftCenter<uint8_t>(tripleddftFinalImageData);
        myImageWriter->write(finalDFT, output_filename, myHeaderImage);
    } else {

        const int N = soundData[0].cols();

        std::shared_ptr<LineOnlyBlueSteinMethod> myBS = std::make_shared<LineOnlyBlueSteinMethod>();
        std::shared_ptr<LineOnlyInvBlueSteinMethod> myIBS = std::make_shared<LineOnlyInvBlueSteinMethod>();

        std::vector<int> factors = myBS->getFactors(N);

        std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> matrixData = myBS->vectorToMatrix<uint16_t>(soundData, factors);
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> soundDataDFT = myBS->computeMethod<uint16_t, double>(matrixData);
        
        displayPreparation(255.0, soundDataDFT);
        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> soundDataDFTUint8 = CastMatrix<double, uint8_t>(soundDataDFT);
        
        BMPHeader myHeaderImage;
        initBMPHeader(myHeaderImage, soundDataDFT[0].cols(), soundDataDFT[0].rows(), 24);
        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> tripleddftFinalImageData(3, soundDataDFTUint8[0]);
        std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> finalDFT = myBandFilter->shiftCenter<uint8_t>(tripleddftFinalImageData);
        myImageWriter->write(finalDFT, output_filename, myHeaderImage);
        
    }  
    
}