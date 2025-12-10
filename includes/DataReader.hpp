
#ifndef DATAREADER
#define DATAREADER 


#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <cassert>
#include "Transforms.hpp"
#include "DataTypes.hpp"
#include <fstream>
#include "FormatFactoring.hpp"


template <typename T>
class AbstractReader 
{

public: 

    AbstractReader(const std::string name, bool applyTransforms, bool verbose) : name_(name), applyTransforms_(applyTransforms), 
    verbose_(verbose) {
        ConverterMap["wav"] = std::make_shared<WAVConverter>();
        ConverterMap["mp3"] = std::make_shared<MP3Converter>();
        ConverterMap["bmp"] = std::make_shared<BMPConverter>();
        ConverterMap["png"] = std::make_shared<PNGConverter>();
    };

    virtual std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> getData(std::string file_path) const = 0;
    virtual int getSize(std::string file_path) const = 0;

    void applyTransform(std::shared_ptr<AbstractTransform<T>> transform, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const;

    void applySeveralTransforms(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const;

    void addTransform(std::shared_ptr<AbstractTransform<T>> transform);
    void removeTransform(int transformIndex);
    void reorganizeTransforms(const std::vector<int>& permutation);

    std::string getPreExtension(std::string file_path) const;

    

    void eraseFile(std::string file_path) const;

    std::string getExtension(std::string file_path) const;
    
    virtual ~AbstractReader(){};

    static std::map<std::string, std::shared_ptr<AbstractConverter>> ConverterMap;

protected: 

    std::string name_;
    bool applyTransforms_;
    bool verbose_;
    std::vector<std::shared_ptr<AbstractTransform<T>>> transformsVector_;

};

class SoundReader : public AbstractReader<uint16_t>
{
public: 

    SoundReader(bool applyTransforms, bool verbose) : AbstractReader<uint16_t>::AbstractReader("Sound Reader", applyTransforms, verbose) {}
    virtual std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> getData(std::string file_path) const override;
    WAVHeader getHeader(std::string file_path) const;
    WAVHeader getFinalHeader(std::string file_path) const;
    ~SoundReader(){};

private: 

    int maxByteLength_;
    std::vector<char> getRawData(std::string file_path) const;
    const int getNumChannels(std::string file_path) const;
    const int getDataSize(std::string file_path) const;
    virtual int getSize(std::string file_path) const override;


};

class ImageReader : public AbstractReader<uint8_t>
{
public: 
    ImageReader(bool applyTransforms, bool verbose) : AbstractReader<uint8_t>::AbstractReader("Image Reader", applyTransforms, verbose) {}
    virtual std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> getData(std::string file_path) const override;
    BMPHeader getHeader(std::string file_path) const;
    BMPHeader getFinalHeader(std::string file_path) const;
    ~ImageReader(){};

private: 

    int maxByteLength_;
    std::vector<char> getRawData(std::string file_path) const;
    
    const int getDataSize(std::string file_path) const;
    virtual int getSize(std::string file_path) const override;
    const int getHeight(std::string file_path) const;
    const int getWidth(std::string file_path) const;
    
};

// Abstract Reader

template <typename T>
std::map<std::string, std::shared_ptr<AbstractConverter>> AbstractReader<T>::ConverterMap;

template <typename T>
std::string AbstractReader<T>::getExtension(std::string file_path) const{
    size_t pos = file_path.find_last_of('.');
    return file_path.substr(pos+1, file_path.size());
}

template <typename T>
void AbstractReader<T>::applyTransform(std::shared_ptr<AbstractTransform<T>> transform, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const {
    transform->apply(matrix);
}

template <typename T>
void AbstractReader<T>::applySeveralTransforms(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const {
    for (auto it = transformsVector_.begin(); it != transformsVector_.end(); it ++){
        (*it)->apply(matrix);
    }
}

template <typename T>
void AbstractReader<T>::addTransform(std::shared_ptr<AbstractTransform<T>> transform) {
    transformsVector_.push_back(transform);
}

template <typename T>
void AbstractReader<T>::removeTransform(int transformIndex){
    assert(transformIndex < transformsVector_.size());
    transformsVector_.erase(transformsVector_.begin() + transformIndex);
}

template <typename T>
void AbstractReader<T>::reorganizeTransforms(const std::vector<int>& permutation){
    assert(permutation.size() == transformsVector_.size());
    std::vector<std::shared_ptr<AbstractTransform<T>>> oldTransVector = transformsVector_;
    for (int i = 0; i < transformsVector_.size(); i ++){
        transformsVector_[i] = oldTransVector[permutation[i]];
    }
}

template <typename T>
void AbstractReader<T>::eraseFile(std::string file_path) const{
    std::string cmd = "rm -rf " + file_path;
    std::system(cmd.c_str());
}

template <typename T>
std::string AbstractReader<T>::getPreExtension(std::string file_path) const {
    size_t pos = file_path.find_last_of('.');
    return file_path.substr(0, pos);
}

///////////////////////////////////////////

#endif