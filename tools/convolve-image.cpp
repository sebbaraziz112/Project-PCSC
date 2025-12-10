#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/methods.hpp"
#include "../includes/utils.hpp"

int main(int argc, char* argv[]){
    if (argc < 3) {
        std::cout << "Usage: imsonpro convolve-image <convolve> <file>\n";
        std::cout << "Type imsonpro to see the details";
        return 0;
    }

    std::map<std::string, std::function<std::shared_ptr<AbstractMethod>()>> factory;
    factory["SobelX"] = [](){ return std::make_shared<SobelX>(); };
    factory["SobelY"] = [](){ return std::make_shared<SobelY>(); };
    factory["Blurr"] = [](){ return std::make_shared<BlurrGaussianConv<10>>(10.0, 10.0); };
    factory["LapClass"] = [](){ return std::make_shared<LaplacianClassical>(); };
    factory["Lap4"] = [](){ return std::make_shared<Laplacian4Connected>(); };
    factory["Lap8"] = [](){ return std::make_shared<Laplacian8Connected>(); };
    factory["lapGauss"] = [](){ return std::make_shared<LaplacianOfGaussian>(); };

    std::shared_ptr<ImageReader> myImageReader = std::make_shared<ImageReader>(true, true);
    std::shared_ptr<ImageWriter> myImageWriter = std::make_shared<ImageWriter>();

    std::string convolve = argv[1];
    std::string input_filename = argv[2];

    std::string output_filename = myImageReader->getPreExtension(input_filename) + "_" + convolve +".bmp";

    std::shared_ptr<AbstractMethod> myConv = factory[convolve]();

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> imageData = myImageReader->getData(input_filename);
    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> imageDataConvolved = myConv->computeMethod<uint8_t, uint8_t>(imageData);


    BMPHeader myHeaderImage = myImageReader->getFinalHeader(input_filename);
    myImageWriter->write(imageDataConvolved, output_filename, myHeaderImage);


}