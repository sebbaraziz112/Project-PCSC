#include <iostream>
#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/methods.hpp"
#include "../includes/utils.hpp"

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
    std::shared_ptr<BlueSteinMethod> myBS = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<InvBlueSteinMethod> myIBS = std::make_shared<InvBlueSteinMethod>();
    std::shared_ptr<BandFiltering> myBandFilter = std::make_shared<BandFiltering>();

    myBandFilter->setLowerFractionX(prct_bottom / 100.0);
    myBandFilter->setUpperFractionX(prct_top / 100.0);
    myBandFilter->setLowerFractionY(prct_bottom / 100.0);
    myBandFilter->setUpperFractionY(prct_top / 100.0);
    myBandFilter->setAmplificationFactor(1.0);

    std::string output_filename = mySoundReader->getPreExtension(input_filename) + "_bandfiltered.wav";

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundData = mySoundReader->getData(input_filename);
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> soundDataDFT = myBS->computeMethod<uint16_t, std::complex<double>>(soundData);
    //std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> bandFilteredDFT = myBandFilter->computeMethod<std::complex<double>, std::complex<double>>(soundDataDFT);
    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundDataIDFT = myIBS->computeMethod<std::complex<double>, uint16_t>(soundDataDFT);
    
    //displayPreparation(65535.0, soundDataIDFT);
    //std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> finalSoundData = CastMatrix<double, uint16_t>(soundDataIDFT);
    WAVHeader myHeaderSound = mySoundReader->getFinalHeader(input_filename);
    mySoundWriter->write(soundData, output_filename, myHeaderSound);

}