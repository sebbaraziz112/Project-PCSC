#ifndef FORMATFACTORING
#define FORMATFACTORING


#include <iostream>
#include <string>
#include <cstdlib>
#include <cassert>
#include <cmath>



class AbstractConverter
{

public: 

    AbstractConverter(std::string name): name_(name){};

    virtual std::string convert(std::string input_file, std::string output_file = "") const = 0;

    std::string createOutput(std::string input_file, std::string sourceExt, std::string destExt, std::string output_file = "") const;

    virtual ~AbstractConverter(){};

private: 

    std::string name_;

};

class WAVConverter: public AbstractConverter
{

public: 

    WAVConverter(): AbstractConverter("WAV to WAV converter"){};

    virtual std::string convert(std::string input_file, std::string output_file = "") const override;

    ~WAVConverter(){};

private: 

    std::string name_;

};

class MP3Converter: public AbstractConverter
{

public: 

    MP3Converter(): AbstractConverter("MP3 to WAV converter"){};

    virtual std::string convert(std::string input_file, std::string output_file = "") const override;

    ~MP3Converter(){};

private: 

    std::string name_;

};

class BMPConverter: public AbstractConverter
{

public: 

    BMPConverter(): AbstractConverter("BMP to BMP converter"){};

    virtual std::string convert(std::string input_file, std::string output_file = "") const override;

    ~BMPConverter(){}

private: 

    std::string name_;

};

class PNGConverter: public AbstractConverter
{

public: 

    PNGConverter(): AbstractConverter("PNG to BMP converter"){};

    virtual std::string convert(std::string input_file, std::string output_file = "") const override;

    ~PNGConverter(){}

private: 

    std::string name_;

};

#endif 



