#include "../includes/FormatFactoring.hpp"


std::string AbstractConverter::createOutput(std::string input_file, std::string sourceExt, std::string destExt, std::string output_file) const{
    size_t inPos = input_file.find_last_of('.');
    assert(input_file.substr(inPos+1, input_file.size()) == sourceExt);
    if (output_file == "") {
        size_t pos = input_file.find_last_of('.');
        output_file = input_file.substr(0, pos) + "." + destExt;
    } else {
        size_t pos = output_file.find_last_of('.');
        assert(output_file.substr(pos+1, output_file.size()) == destExt);
    }
    return output_file;

}

std::string WAVConverter::convert(std::string input_file, std::string output_file) const {
    return input_file;
}

std::string MP3Converter::convert(std::string input_file, std::string output_file) const {
    output_file = AbstractConverter::createOutput(input_file, "mp3", "wav", output_file);
    std::string cmd = "mpg123 -w " + output_file + " " + input_file +  " > /dev/null 2>&1";
    std::system(cmd.c_str());
    return output_file;
}

std::string BMPConverter::convert(std::string input_file, std::string output_file) const {
    output_file = AbstractConverter::createOutput(input_file, "bmp", "bmp", output_file);
    std::string cmd = "convert " + input_file +  " -type TrueColor " +  "bmp3:" + output_file;
    std::system(cmd.c_str());
    return output_file;
}

std::string PNGConverter::convert(std::string input_file, std::string output_file) const {
    output_file = AbstractConverter::createOutput(input_file, "png", "bmp", output_file);
    std::string cmd = "convert " + input_file +  " -type TrueColor " +  "bmp3:" + output_file;
    std::system(cmd.c_str());
    return output_file;
}