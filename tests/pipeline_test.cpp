#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/methods.hpp"
#include "../includes/utils.hpp"


int main(int argc, char* argv[]){

    std::string original_file_image = "../ressources/ImageFolder/sine_diag.png";
    std::string new_file_image = "../ressources/ImageFolder/test.bmp";

    std::shared_ptr<ImageReader> myImageReader = std::make_shared<ImageReader>(true, true);
    std::shared_ptr<ImageWriter> myImageWriter = std::make_shared<ImageWriter>();
    
    std::shared_ptr<BlueSteinMethod> myBS = std::make_shared<BlueSteinMethod>();
    std::shared_ptr<InvBlueSteinMethod> myIBS = std::make_shared<InvBlueSteinMethod>();

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> imageData = myImageReader->getData(original_file_image);
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> imageDataDFT = myBS->computeMethod<uint8_t, double>(imageData);
    centerAround(255.0, imageDataDFT);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> newData = CastMatrix<double, uint8_t>(imageDataDFT);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> shiftedData = shiftCenter<uint8_t>(newData);

    BMPHeader myHeaderImage = myImageReader->getFinalHeader(original_file_image);
    myImageWriter->write(shiftedData, new_file_image, myHeaderImage);


}



