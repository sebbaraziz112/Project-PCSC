#include <iostream>
#include "../includes/Visual.hpp"
#include "../includes/methods.hpp"
#include "../includes/DataReader.hpp"

int main(int argc, char* argv[]) {
    
    if (argc < 2) {
        std::cout << "Usage: imsonpro histogram-sound <file>\n";
        return 0;
    }

    std::string filePath = argv[1];
    

    std::shared_ptr<SoundReader> mySoundReader = std::make_shared<SoundReader>(true, true);
    std::string output_filename = mySoundReader->getPreExtension(filePath) + "_histogram";
    

    std::shared_ptr<ProbDensity> probDensity = std::make_shared<ProbDensity>();
    std::shared_ptr<Histogram> histogramPlot = std::make_shared<Histogram>();

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundData = mySoundReader->getData(filePath);
    std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> soundDataPDF = probDensity->computeMethod<uint16_t, int>(soundData);

    for (int channel = 0; channel < soundDataPDF.size(); channel ++){
        histogramPlot->setTitle("Channel " + std::to_string(channel) + " Histogram");
        histogramPlot->setAxisXLabel("Intensity");
        histogramPlot->setAxisYLabel("Frequency");
        histogramPlot->plot(soundDataPDF[channel], output_filename + "_channel_" + std::to_string(channel) + ".png");
    }
 
}