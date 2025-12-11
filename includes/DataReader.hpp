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

/**
 * @brief Abstract class defining a generic data reader.
 *
 * This class handles extensions, format conversion, transform application,
 * and supports multiple file types (audio, image).
 *
 * @tparam T Scalar type of the data read (e.g., uint8_t for images, uint16_t for audio).
 */
template <typename T>
class AbstractReader 
{

public: 

    /**
     * @brief Main constructor.
     *
     * Initializes the reader name, transform application flag, verbosity,
     * and registers the supported file converters.
     *
     * @param name Name of the reader.
     * @param applyTransforms Whether transforms should be applied automatically.
     * @param verbose Enables or disables verbose mode.
     */
    AbstractReader(const std::string name, bool applyTransforms, bool verbose) : name_(name), applyTransforms_(applyTransforms), 
    verbose_(verbose) {
        ConverterMap["wav"] = std::make_shared<WAVConverter>();
        ConverterMap["mp3"] = std::make_shared<MP3Converter>();
        ConverterMap["bmp"] = std::make_shared<BMPConverter>();
        ConverterMap["png"] = std::make_shared<PNGConverter>();
    };

    /**
     * @brief Reads a file and returns its content as a vector of Eigen matrices.
     *
     * @param file_path Path to the file.
     * @return Vector of matrices containing the extracted data.
     */
    virtual std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> getData(std::string file_path) const = 0;

    /**
     * @brief Returns the raw size of the file in bytes.
     *
     * @param file_path Path to the file.
     * @return File size in bytes.
     */
    virtual int getSize(std::string file_path) const = 0;

    /**
     * @brief Applies a single transform to a matrix.
     *
     * @param transform Transform to apply.
     * @param matrix Matrix to transform.
     */
    void applyTransform(std::shared_ptr<AbstractTransform<T>> transform, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const;

    /**
     * @brief Applies all registered transforms to a matrix.
     *
     * @param matrix Matrix to transform.
     */
    void applySeveralTransforms(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &matrix) const;

    /**
     * @brief Adds a transform to the internal list.
     *
     * @param transform Transform to add.
     */
    void addTransform(std::shared_ptr<AbstractTransform<T>> transform);

    /**
     * @brief Removes a transform by index.
     *
     * @param transformIndex Index in the transform vector.
     */
    void removeTransform(int transformIndex);

    /**
     * @brief Reorders transforms based on a given permutation.
     *
     * @param permutation Vector specifying the new order.
     */
    void reorganizeTransforms(const std::vector<int>& permutation);

    /**
     * @brief Returns the file path without its extension.
     *
     * @param file_path Input file path.
     * @return File path without extension.
     */
    std::string getPreExtension(std::string file_path) const;

    /**
     * @brief Deletes a file on disk.
     *
     * @param file_path Path to the file to remove.
     */
    void eraseFile(std::string file_path) const;

    /**
     * @brief Extracts the extension from a file name.
     *
     * @param file_path Path to the file.
     * @return Extracted extension.
     */
    std::string getExtension(std::string file_path) const;

    virtual ~AbstractReader(){};

    /**
     * @brief Static map linking file extensions to their converters.
     */
    static std::map<std::string, std::shared_ptr<AbstractConverter>> ConverterMap;

protected: 

    /** @brief Name of the reader (e.g., "Sound Reader"). */
    std::string name_;

    /** @brief Whether transforms should be applied when reading. */
    bool applyTransforms_;

    /** @brief Enable or disable verbose output. */
    bool verbose_;

    /** @brief List of transforms applied to the data. */
    std::vector<std::shared_ptr<AbstractTransform<T>>> transformsVector_;

};

/**
 * @brief Reader specialized for audio files (e.g., WAV, MP3).
 */
class SoundReader : public AbstractReader<uint16_t>
{
public: 

    /**
     * @brief Constructor for the audio reader.
     *
     * @param applyTransforms Whether to automatically apply transforms.
     * @param verbose Whether verbose mode is enabled.
     */
    SoundReader(bool applyTransforms, bool verbose) : AbstractReader<uint16_t>::AbstractReader("Sound Reader", applyTransforms, verbose) {}

    /**
     * @brief Reads an audio file and returns the data matrices.
     *
     * @param file_path Path to the audio file.
     */
    virtual std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> getData(std::string file_path) const override;

    /**
     * @brief Retrieves the raw WAV header.
     *
     * @param file_path Path to the WAV file.
     */
    WAVHeader getHeader(std::string file_path) const;

    /**
     * @brief Retrieves the final WAV header after conversion.
     *
     * @param file_path Path to the WAV file.
     */
    WAVHeader getFinalHeader(std::string file_path) const;

    ~SoundReader(){};

private: 

    /** @brief Maximum number of bytes to read from the file. */
    int maxByteLength_;

    /**
     * @brief Reads the raw audio byte content.
     *
     * @param file_path Path to the file.
     * @return Raw data buffer.
     */
    std::vector<char> getRawData(std::string file_path) const;

    /**
     * @brief Returns the number of audio channels.
     */
    const int getNumChannels(std::string file_path) const;

    /**
     * @brief Returns the size of the audio data section.
     */
    const int getDataSize(std::string file_path) const;

    /**
     * @brief Returns the raw file size in bytes.
     */
    virtual int getSize(std::string file_path) const override;

};

/**
 * @brief Reader specialized for image files (e.g., BMP, PNG).
 */
class ImageReader : public AbstractReader<uint8_t>
{
public: 

    /**
     * @brief Constructor for the image reader.
     */
    ImageReader(bool applyTransforms, bool verbose) : AbstractReader<uint8_t>::AbstractReader("Image Reader", applyTransforms, verbose) {}

    /**
     * @brief Reads an image file and returns matrices representing its data.
     *
     * @param file_path Path to the image file.
     */
    virtual std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> getData(std::string file_path) const override;

    /**
     * @brief Retrieves the raw BMP header.
     *
     * @param file_path Path to the file.
     */
    BMPHeader getHeader(std::string file_path) const;

    /**
     * @brief Retrieves the final BMP header after processing.
     *
     * @param file_path Path to the file.
     */
    BMPHeader getFinalHeader(std::string file_path) const;

    ~ImageReader(){};

private: 

    /** @brief Maximum data size in bytes. */
    int maxByteLength_;

    /**
     * @brief Returns the raw pixel data of the image.
     */
    std::vector<char> getRawData(std::string file_path) const;
    
    /**
     * @brief Returns the size of the image data block.
     */
    const int getDataSize(std::string file_path) const;

    /**
     * @brief Returns the raw file size in bytes.
     */
    virtual int getSize(std::string file_path) const override;

    /**
     * @brief Returns the image height.
     */
    const int getHeight(std::string file_path) const;

    /**
     * @brief Returns the image width.
     */
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
