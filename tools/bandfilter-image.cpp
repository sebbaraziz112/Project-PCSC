#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/methods.hpp"
#include "../includes/utils.hpp"



int main(int argc, char* argv[]){
    if (argc < 4) {
        std::cout << "Usage: imsonpro bandfilter-image <file> <prct bottom> <prct top>\n";
        return 0;
    }

    std::string input_filename = argv[1];
    double prct_bottom = std::stod(argv[2]);
    double prct_top = std::stod(argv[3]);


    std::shared_ptr<ImageReader> myImageReader = std::make_shared<ImageReader>(true, true);
    std::shared_ptr<ImageWriter> myImageWriter = std::make_shared<ImageWriter>();
    std::shared_ptr<BlueSteinMethod> myBS = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<InvBlueSteinMethod> myIBS = std::make_shared<InvBlueSteinMethod>();
    std::shared_ptr<BandFiltering> myBandFilter = std::make_shared<BandFiltering>();

    myBandFilter->setLowerFractionX(prct_bottom / 100.0);
    myBandFilter->setUpperFractionX(prct_top / 100.0);
    myBandFilter->setLowerFractionY(prct_bottom / 100.0);
    myBandFilter->setUpperFractionY(prct_top / 100.0);
    myBandFilter->setAmplificationFactor(1.0);

    std::string output_filename = myImageReader->getPreExtension(input_filename) + "_bandfiltered.bmp";
    std::string bf_dft_filename = myImageReader->getPreExtension(input_filename) + "_bandfiltered_dft.bmp";
    std::string dft_filename = myImageReader->getPreExtension(input_filename) + "_dft.bmp";

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> imageData = myImageReader->getData(input_filename);
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> imageDataDFT = myBS->computeMethod<uint8_t, std::complex<double>>(imageData);
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> shiftedDataDFT = myBandFilter->shiftCenter<std::complex<double>>(imageDataDFT);
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> bandFilteredDFT = myBandFilter->computeMethod<std::complex<double>, std::complex<double>>(shiftedDataDFT);
    std::vector<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> invShiftedDataDFT = myBandFilter->shiftInverse<std::complex<double>>(bandFilteredDFT);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> imageDataIDFT = myIBS->computeMethod<std::complex<double>, uint8_t>(invShiftedDataDFT);
    

    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> bfDftData = CastMatrix<std::complex<double>, double>(bandFilteredDFT);
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> dftData = CastMatrix<std::complex<double>, double>(imageDataDFT);
    displayPreparation(255.0, bfDftData);
    displayPreparation(255.0, dftData);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> bfDftFinalData = CastMatrix<double, uint8_t>(bfDftData);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> dftFinalData = CastMatrix<double, uint8_t>(dftData);


    BMPHeader myHeaderImage = myImageReader->getFinalHeader(input_filename);
    myImageWriter->write(imageDataIDFT, output_filename, myHeaderImage);
    myImageWriter->write(bfDftFinalData, bf_dft_filename, myHeaderImage);
    myImageWriter->write(dftFinalData, dft_filename, myHeaderImage);

}