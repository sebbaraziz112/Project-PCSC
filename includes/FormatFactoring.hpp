/**
 * @file FormatFactoring.hpp
 * @brief Defines abstract and concrete classes for file format conversion (audio and image).
 *
 * This header provides the AbstractConverter base class, which defines a generic
 * interface for file conversion, and several concrete converters for WAV, MP3, BMP, and PNG formats.
 */

#ifndef FORMATFACTORING
#define FORMATFACTORING

#include <iostream>
#include <string>
#include <cstdlib>
#include <cassert>
#include <cmath>

/**
 * @class AbstractConverter
 * @brief Abstract base class for file format converters.
 *
 * Provides a generic interface for converting input files to a specific format.
 */
class AbstractConverter
{
public:

    /**
     * @brief Constructor.
     * @param name Name of the converter.
     */
    AbstractConverter(std::string name) : name_(name) {};

    /**
     * @brief Convert a file from one format to another.
     * @param input_file Path to the input file.
     * @param output_file Optional path to the output file.
     * @return Path of the converted file as a string.
     *
     * This is a pure virtual function that must be implemented in derived classes.
     */
    virtual std::string convert(std::string input_file, std::string output_file = "") const = 0;

    /**
     * @brief Helper function to create an output file name based on input, source, and destination extensions.
     * @param input_file Path to the input file.
     * @param sourceExt Source file extension.
     * @param destExt Destination file extension.
     * @param output_file Optional output file path.
     * @return Path of the output file as a string.
     */
    std::string createOutput(std::string input_file, std::string sourceExt, std::string destExt, std::string output_file = "") const;

    /**
     * @brief Destructor.
     */
    virtual ~AbstractConverter() {};

private:
    std::string name_; /**< Name of the converter */
};

/**
 * @class WAVConverter
 * @brief Converts WAV files to WAV format (e.g., for preprocessing or normalization).
 */
class WAVConverter : public AbstractConverter
{
public:
    WAVConverter() : AbstractConverter("WAV to WAV converter") {};

    virtual std::string convert(std::string input_file, std::string output_file = "") const override;

    ~WAVConverter() {};

private:
    std::string name_;
};

/**
 * @class MP3Converter
 * @brief Converts MP3 files to WAV format.
 */
class MP3Converter : public AbstractConverter
{
public:
    MP3Converter() : AbstractConverter("MP3 to WAV converter") {};

    virtual std::string convert(std::string input_file, std::string output_file = "") const override;

    ~MP3Converter() {};

private:
    std::string name_;
};

/**
 * @class BMPConverter
 * @brief Converts BMP files to BMP format (e.g., for preprocessing or normalization).
 */
class BMPConverter : public AbstractConverter
{
public:
    BMPConverter() : AbstractConverter("BMP to BMP converter") {};

    virtual std::string convert(std::string input_file, std::string output_file = "") const override;

    ~BMPConverter() {};

private:
    std::string name_;
};

/**
 * @class PNGConverter
 * @brief Converts PNG files to BMP format.
 */
class PNGConverter : public AbstractConverter
{
public:
    PNGConverter() : AbstractConverter("PNG to BMP converter") {};

    virtual std::string convert(std::string input_file, std::string output_file = "") const override;

    ~PNGConverter() {};

private:
    std::string name_;
};

#endif // FORMATFACTORING
