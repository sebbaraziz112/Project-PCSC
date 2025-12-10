#ifndef DATAWRITER
#define DATAWRITER


#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <cassert>
#include <fstream>
#include <cmath>
#include "DataTypes.hpp"
#include "DataReader.hpp"



template <typename T, typename H>
class AbstractWriter
{

public:

    AbstractWriter(std::string name, bool verbose): name_(name), verbose_(verbose){};
    virtual ~AbstractWriter(){};

    virtual void writeData(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& data, H& header, std::fstream& file_output) const = 0;

    void writeHeader(H& header, std::fstream& file_output) const;
    std::string write(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& data, std::string filename, H header) const;
    std::string getPreExtension(std::string file_path) const;
    std::string getExtension(std::string file_path) const;

    

protected:

    std::string name_;
    bool verbose_;
    std::string extension_;


};

class SoundWriter: public AbstractWriter<uint16_t, WAVHeader>
{

public: 

    SoundWriter(): AbstractWriter<uint16_t, WAVHeader>("Sound Writer", true){
        extension_ = "wav";
    };
    ~SoundWriter(){};
    void writeData(std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>>& data, WAVHeader& header, std::fstream& file_output) const override;

private:

};

class ImageWriter: public AbstractWriter<uint8_t, BMPHeader>
{

public: 

    ImageWriter(): AbstractWriter<uint8_t, BMPHeader>("Image Writer", true){
        extension_ = "bmp";
    };
    ~ImageWriter(){};
    void writeData(std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>>& data, BMPHeader& header, std::fstream& file_output) const override;

private:

};


// Abstract Writer

template <typename T, typename H>
void AbstractWriter<T, H>::writeHeader(H& header, std::fstream& file_output) const{
    file_output.write(reinterpret_cast<const char*>(&header), sizeof(header));
}

template <typename T, typename H>
std::string  AbstractWriter<T, H>::write(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& data, std::string filename, H header) const{
    std::fstream file_output(this->getPreExtension(filename) + "." + extension_, std::ios::out | std::ios::trunc);
    this->writeHeader(header, file_output);
    this->writeData(data, header, file_output);
    file_output.close();
    return this->getPreExtension(filename) + "." + extension_;
}

template <typename T, typename H>
std::string AbstractWriter<T, H>::getExtension(std::string file_path) const {
    size_t pos = file_path.find_last_of('.');
    return file_path.substr(pos+1, file_path.size());
}

template <typename T, typename H>
std::string AbstractWriter<T, H>::getPreExtension(std::string file_path) const {
    size_t pos = file_path.find_last_of('.');
    return file_path.substr(0, pos);
}

/////////////////////////////////////////////


#endif