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
    

    std::shared_ptr<ImageReader> myImageReader = std::make_shared<ImageReader>(true, true);
    std::string output_filename_r = myImageReader->getPreExtension(filePath) + "_histogram_R";
    std::string output_filename_g = myImageReader->getPreExtension(filePath) + "_histogram_G";
    std::string output_filename_b = myImageReader->getPreExtension(filePath) + "_histogram_B";

    std::shared_ptr<ProbDensity> probDensity = std::make_shared<ProbDensity>();
    std::shared_ptr<Histogram> histogramPlot = std::make_shared<Histogram>();

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> imageData = myImageReader->getData(filePath);
    std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> imageDataPDF = probDensity->computeMethod<uint8_t, int>(imageData);

    histogramPlot->setTitle("Red Channel Histogram");
    histogramPlot->setAxisXLabel("Intensity");
    histogramPlot->setAxisYLabel("Frequency");
    histogramPlot->plot(imageDataPDF[0], output_filename_r + ".png");

    histogramPlot->setTitle("Green Channel Histogram");
    histogramPlot->setAxisXLabel("Intensity");
    histogramPlot->setAxisYLabel("Frequency");
    histogramPlot->plot(imageDataPDF[1], output_filename_g + ".png");

    histogramPlot->setTitle("Blue Channel Histogram");
    histogramPlot->setAxisXLabel("Intensity");
    histogramPlot->setAxisYLabel("Frequency");
    histogramPlot->plot(imageDataPDF[2], output_filename_b + ".png");
 
}