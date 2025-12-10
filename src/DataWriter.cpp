#include "../includes/DataWriter.hpp"
#include "../includes/DataTypes.hpp"


// Abstract Writer

/////////////////////////////////////////////

//  Sound Writer

void SoundWriter::writeData(std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>>& data, WAVHeader& header, std::fstream& file_output) const {
    const int nCols = data[0].cols();
    const int nChann = data.size();
    std::vector<uint16_t> reformatData;
    for (int i = 0; i < nCols; i ++){
        for (int chann = 0; chann < nChann; chann ++){
            reformatData.push_back(data[chann](0, i));
        }
    }
    file_output.write(reinterpret_cast<char*>(reformatData.data()), reformatData.size() * sizeof(uint16_t));
}

/////////////////////////////////////////////

// Image Writer

void ImageWriter::writeData(std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>>& data, BMPHeader& header, std::fstream& file_output) const {

    const int nRows = data[0].rows();
    const int nCols = data[0].cols();
    const int nChann = 3;

    const int paddingBytes = (4 - ((nCols * 3) % 4)) % 4;

    std::vector<uint8_t> reformatData;

    for (int i = 0; i < nRows; i ++){
        for (int j = 0; j < nCols; j ++){
            for (int chann = 0; chann < nChann; chann ++){
                reformatData.push_back(data[chann](i, j));
            }
        }
        for (int k = 0; k < paddingBytes; k ++){
            reformatData.push_back(static_cast<uint8_t>(0));
        }
    }
    file_output.write(reinterpret_cast<char*>(reformatData.data()), reformatData.size() * sizeof(uint8_t));
}

/////////////////////////////////////////////