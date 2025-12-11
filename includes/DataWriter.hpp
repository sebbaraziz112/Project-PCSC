/**
 * @file DataWriter.hpp
 * @brief Defines abstract and concrete classes for writing audio and image data to files.
 *
 * Provides the AbstractWriter template for generic data writing and concrete implementations
 * for WAV (SoundWriter) and BMP (ImageWriter) formats.
 */

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

/**
 * @class AbstractWriter
 * @brief Abstract base class template for writing matrices to files.
 *
 * @tparam T The type of matrix elements (e.g., uint8_t, uint16_t).
 * @tparam H The type of the header (e.g., WAVHeader, BMPHeader).
 */
template <typename T, typename H>
class AbstractWriter
{
public:
    /**
     * @brief Constructor.
     * @param name Name of the writer.
     * @param verbose If true, prints progress messages.
     */
    AbstractWriter(std::string name, bool verbose) : name_(name), verbose_(verbose) {};

    /**
     * @brief Virtual destructor.
     */
    virtual ~AbstractWriter() {};

    /**
     * @brief Pure virtual method to write data matrices to a file.
     * @param data Vector of matrices to write.
     * @param header Header containing metadata.
     * @param file_output Opened file stream for writing.
     */
    virtual void writeData(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& data, H& header, std::fstream& file_output) const = 0;

    /**
     * @brief Writes the header to the output file.
     * @param header Header to write.
     * @param file_output Opened file stream.
     */
    void writeHeader(H& header, std::fstream& file_output) const;

    /**
     * @brief Writes data and header to a file.
     * @param data Vector of matrices to write.
     * @param filename Output file path.
     * @param header Header to write.
     * @return The actual output file path used.
     */
    std::string write(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& data, std::string filename, H header) const;

    /**
     * @brief Extracts the extension from a file path.
     * @param file_path File path string.
     * @return File extension as string.
     */
    std::string getExtension(std::string file_path) const;

    /**
     * @brief Extracts the file path without extension.
     * @param file_path File path string.
     * @return File path excluding extension.
     */
    std::string getPreExtension(std::string file_path) const;

protected:
    std::string name_;   /**< Name of the writer */
    bool verbose_;       /**< Verbose flag */
    std::string extension_; /**< File extension used for output */
};

/**
 * @class SoundWriter
 * @brief Concrete writer for WAV files (audio).
 */
class SoundWriter : public AbstractWriter<uint16_t, WAVHeader>
{
public:
    SoundWriter() : AbstractWriter<uint16_t, WAVHeader>("Sound Writer", true) { extension_ = "wav"; };
    ~SoundWriter() {};

    /**
     * @brief Write audio data matrices to WAV file.
     */
    void writeData(std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>>& data, WAVHeader& header, std::fstream& file_output) const override;

private:
};

/**
 * @class ImageWriter
 * @brief Concrete writer for BMP files (images).
 */
class ImageWriter : public AbstractWriter<uint8_t, BMPHeader>
{
public:
    ImageWriter() : AbstractWriter<uint8_t, BMPHeader>("Image Writer", true) { extension_ = "bmp"; };
    ~ImageWriter() {};

    /**
     * @brief Write image data matrices to BMP file.
     */
    void writeData(std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>>& data, BMPHeader& header, std::fstream& file_output) const override;

private:
};

// ---------------- Implementation of template methods ----------------

template <typename T, typename H>
void AbstractWriter<T, H>::writeHeader(H& header, std::fstream& file_output) const {
    file_output.write(reinterpret_cast<const char*>(&header), sizeof(header));
}

template <typename T, typename H>
std::string AbstractWriter<T, H>::write(std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& data, std::string filename, H header) const {
    std::fstream file_output(this->getPreExtension(filename) + "." + extension_, std::ios::out | std::ios::trunc);
    this->writeHeader(header, file_output);
    this->writeData(data, header, file_output);
    file_output.close();
    return this->getPreExtension(filename) + "." + extension_;
}

template <typename T, typename H>
std::string AbstractWriter<T, H>::getExtension(std::string file_path) const {
    size_t pos = file_path.find_last_of('.');
    return file_path.substr(pos + 1, file_path.size());
}

template <typename T, typename H>
std::string AbstractWriter<T, H>::getPreExtension(std::string file_path) const {
    size_t pos = file_path.find_last_of('.');
    return file_path.substr(0, pos);
}

#endif // DATAWRITER
