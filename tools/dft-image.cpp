#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/methods.hpp"
#include "../includes/utils.hpp"

int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Usage: imsonpro dft-image <file>\n";
        return 0;
    }

    std::string input_filename = argv[1];
    std::shared_ptr<ImageReader> myImageReader = std::make_shared<ImageReader>(true, true);
    std::shared_ptr<ImageWriter> myImageWriter = std::make_shared<ImageWriter>();
    std::shared_ptr<BlueSteinMethod> myBS = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<BandFiltering> myBandFilter = std::make_shared<BandFiltering>();

    std::string output_filename = myImageReader->getPreExtension(input_filename) + "_dft.bmp";


    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> imageData = myImageReader->getData(input_filename);
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> imageDataDFT = myBS->computeMethod<uint8_t, double>(imageData);
    displayPreparation(255.0, imageDataDFT);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> finalData = CastMatrix<double, uint8_t>(imageDataDFT);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> shiftedFinalData = myBandFilter->shiftCenter<uint8_t>(finalData);

    BMPHeader myHeaderImage = myImageReader->getFinalHeader(input_filename);
    myImageWriter->write(shiftedFinalData, output_filename, myHeaderImage);


}