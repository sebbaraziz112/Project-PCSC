#include "includes/DataReader.hpp"
#include <iomanip>
#include "includes/methods.hpp"
#include "includes/DataReader.hpp"
#include "includes/DataWriter.hpp"



int main(int argc, char *argv[]){

    if (argc < 2) {
        std::cout << "Usage: imsonpro <command> [args]" << std::endl;
        std::cout << "The availabe commands are: " << std::endl;
        std::cout << "      - dft-image [filename]" << std::endl;
        std::cout << "      - convolve-image [convolve type] [filename]" << std::endl;
        std::cout << "            convolve: -SobelX" << std::endl;
        std::cout << "                      -SobelY" << std::endl;
        std::cout << "                      -Blurr" << std::endl;
        std::cout << "                      -LapClass" << std::endl;
        std::cout << "                      -Lap4" << std::endl;
        std::cout << "                      -Lap8" << std::endl;
        std::cout << "                      -LapGauss" << std::endl;
        std::cout << "      - bandfilter-image [filename] [prct bottom] [prct ceil]" << std::endl;
        std::cout << "      - dft-sound [fileame]" << std::endl;
        std::cout << "      - bandfilter-sound [filename] [prct bottom] [prct ceil]" << std::endl;
        return 0;
    }

    std::string cmd = argv[1];

    std::string executable = "/usr/local/bin/imsonpro_build/imsonpro-" + cmd;
    std::string fullCmd = executable;
    for (int i = 2; i < argc; ++i) {
        fullCmd += " ";
        fullCmd += argv[i];
    }

    int ret = std::system(fullCmd.c_str());
    return ret;

}