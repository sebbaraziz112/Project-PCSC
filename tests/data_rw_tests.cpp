#include "../includes/DataReader.hpp"
#include "../includes/DataWriter.hpp"
#include "../includes/utils.hpp"
#include <gtest/gtest.h>
#include <fstream>


TEST(DATA_READ_WRITE, FULL_PIPELINE_SOUND){

    std::shared_ptr<SoundReader> mySoundReader = std::make_shared<SoundReader>(true, true);
    std::shared_ptr<SoundWriter> mySoundWriter = std::make_shared<SoundWriter>();

    std::string original_file_sound = "../ressources/TrackFolder/track1.mp3";
    std::string new_file_sound = "../ressources/TrackFolder/new_track1.wav";
    std::string generated_file_sound = "../ressources/TrackFolder/track1.wav";


    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> soundData = mySoundReader->getData(original_file_sound);
    WAVHeader myHeaderSound = mySoundReader->getFinalHeader(original_file_sound);
    mySoundWriter->write(soundData, new_file_sound, myHeaderSound);

    EXPECT_EQ(compareFiles(new_file_sound, generated_file_sound), 1);

}

TEST(DATA_READ_WRITE, FULL_PIPELINE_IMAGE){


    std::shared_ptr<ImageReader> myImageReader = std::make_shared<ImageReader>(true, true);
    std::shared_ptr<ImageWriter> myImageWriter = std::make_shared<ImageWriter>();

    std::string original_file_image = "../ressources/ImageFolder/Mauritius_beach.png";
    std::string new_file_image = "../ressources/ImageFolder/new_Mauritius_beach.bmp";
    std::string generated_file_image = "../ressources/ImageFolder/Mauritius_beach.bmp";


    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> imageData = myImageReader->getData(original_file_image);
    BMPHeader myHeaderImage = myImageReader->getFinalHeader(original_file_image);
    myImageWriter->write(imageData, new_file_image, myHeaderImage);

    EXPECT_EQ(compareFiles(new_file_image, generated_file_image), 1);

    

}


