#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/methods.hpp"
#include "../includes/utils.hpp"


int main(int argc, char* argv[]){
    if (argc < 2) {
        std::cout << "Usage: imsonpro dft-sound <file>\n";
        return 0;
    }

    std::string input_filename = argv[1];
    std::shared_ptr<SoundReader> mySoundReader = std::make_shared<SoundReader>(true, true);
    std::shared_ptr<BlueSteinMethod> myBS = std::make_shared<BlueSteinMethod>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundData = mySoundReader->getData(input_filename);
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> soundDataDFT = myBS->computeMethod<uint16_t, double>(soundData);

    std::cout << "Size of the DFT result: " << soundDataDFT[0].rows() << " x " << soundDataDFT[0].cols() << "\n";   
    

}