#include "../includes/DataReader.hpp"

// Sound Reader

WAVHeader SoundReader::getHeader(std::string file_path) const{
    WAVHeader MyHeader;
    std::ifstream file_input(file_path, std::ios::binary);
    file_input.seekg(0, std::ios::beg);
    file_input.read(reinterpret_cast<char*>(&MyHeader), sizeof(MyHeader));
    file_input.close();
    return MyHeader;
}

WAVHeader SoundReader::getFinalHeader(std::string file_path) const {
    WAVHeader MyHeader;
    std::string new_path = (this->ConverterMap[this->getExtension(file_path)])->convert(file_path);
    std::ifstream file_input(new_path, std::ios::binary);
    file_input.seekg(0, std::ios::beg);
    file_input.read(reinterpret_cast<char*>(&MyHeader), sizeof(MyHeader));
    file_input.close();
    return MyHeader;
}

int SoundReader::getSize(std::string file_path) const{
    return this->getHeader(file_path).DataSize;
}

std::vector<char> SoundReader::getRawData(std::string file_path) const{
    std::ifstream file_input(file_path, std::ios::binary);
    file_input.seekg(44, std::ios::beg);
    int DataSize = this->getSize(file_path);
    std::vector<char> vectorChar;
    vectorChar.resize(DataSize);
    file_input.read(vectorChar.data(), DataSize);
    file_input.close();
    return vectorChar;
}

const int SoundReader::getNumChannels(std::string file_path) const {
    return this->getHeader(file_path).NumChannels;
}

const int SoundReader::getDataSize(std::string file_path) const {
    return this->getHeader(file_path).DataSize;
}

std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> SoundReader::getData(std::string file_path) const{
    std::string new_path = (this->ConverterMap[this->getExtension(file_path)])->convert(file_path);
    std::vector<char> MyVector = this->getRawData(new_path);
    uint16_t * RawIntData = reinterpret_cast<uint16_t*>(MyVector.data());
    const int numChan = this->getNumChannels(new_path);
    const int numData = MyVector.size() / sizeof(uint16_t);
    const int numCol = numData / numChan;

    std::vector<Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic>> DataVector;
    for (int i = 0; i < numChan; i ++){
        Eigen::Matrix<uint16_t, Eigen::Dynamic, Eigen::Dynamic> mat(1, numCol);
        DataVector.push_back(mat);
    }

    for (int i = 0; i < numCol; i ++){
        for (int chann = 0; chann < numChan; chann ++){
            DataVector[chann](0, i) = *(RawIntData + (numChan * i + chann));
        }
    }


    if (this->applyTransforms_) {
        for (int i = 0; i < numChan; i ++)
        this->applySeveralTransforms(DataVector[i]);
    }

    if (file_path != new_path){
        this->eraseFile(new_path);
    }

    return DataVector;

}

///////////////////////////////////////////

// Image Reader

BMPHeader ImageReader::getHeader(std::string file_path) const{
    BMPHeader myHeader;
    std::ifstream file_input(file_path, std::ios::binary);
    file_input.seekg(0, std::ios::beg);
    file_input.read(reinterpret_cast<char*>(&myHeader), sizeof(myHeader));
    file_input.close();
    return myHeader;
}

BMPHeader ImageReader::getFinalHeader(std::string file_path) const {
    BMPHeader myHeader;
    std::string new_path = (this->ConverterMap[this->getExtension(file_path)])->convert(file_path);
    std::ifstream file_input(new_path, std::ios::binary);
    file_input.seekg(0, std::ios::beg);
    file_input.read(reinterpret_cast<char*>(&myHeader), sizeof(myHeader));
    file_input.close();
    return myHeader;
}


int ImageReader::getSize(std::string file_path) const {
    return this->getHeader(file_path).FileSize;
}

const int ImageReader::getDataSize(std::string file_path) const{
    return this->getHeader(file_path).ImageSize;
}

const int ImageReader::getHeight(std::string file_path) const {
    return this->getHeader(file_path).Height;
}

const int ImageReader::getWidth(std::string file_path) const {
    return this->getHeader(file_path).Width;
}

std::vector<char> ImageReader::getRawData(std::string file_path) const{
    BMPHeader myHeader = this->getHeader(file_path);
    uint16_t BitsPerPixel = myHeader.BitsPerPixel;
    uint32_t DataOffset = myHeader.DataOffset;
    uint32_t ImageSize = myHeader.ImageSize;
    if (BitsPerPixel != 24){
        std::cerr << "The format is not standard, something is bad !!" << std::endl;
    }

    std::ifstream file_input(file_path, std::ios::binary);
    file_input.seekg(DataOffset, std::ios::beg);

    std::vector<char> data;
    data.resize(ImageSize);
    file_input.read(data.data(), ImageSize);
    return data;
}

std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> ImageReader::getData(std::string file_path) const {
    std::string new_path = (this->ConverterMap[this->getExtension(file_path)])->convert(file_path);
    
    std::vector<char> MyVector = this->getRawData(new_path);

    int Width = this->getWidth(new_path);
    int Height = this->getHeight(new_path);


    int paddingByte = (4 - ((Width * 3) % 4)) % 4;
    int Rowsize = Width * 3 + paddingByte;
    
    int ColumnSize = Height;

    assert(MyVector.size() == Rowsize * ColumnSize);

    // We need to remove some stuff from the MyVector
    int count = 0;
    for (int i = 0; i < Height; i ++){
        int startIdx = Rowsize * i;
        for (int k = 0; k < paddingByte; k ++){
            MyVector.erase(MyVector.begin() + startIdx + Width * 3 + k - count);
            count ++;
        }  
    }

    assert(MyVector.size() % 3 == 0);
    assert(MyVector.size() == Height * Width * 3);

    uint8_t * RawIntData = reinterpret_cast<uint8_t*>(MyVector.data());

    std::vector<Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>> DataVector;

    for (int i = 0; i < 3; i ++){
        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> mat(Height, Width);
        DataVector.push_back(mat);
    }

    for (int i = 0; i < Height; i ++){
        for (int j = 0; j < Width; j ++){

            for (int chann = 0; chann < 3; chann ++){
                DataVector[chann](i, j) = *(RawIntData + (i * Width + j) * 3 + chann);
            }

        }
    }
    

    if (this->applyTransforms_) {
        for (int i = 0; i < 3; i ++){
            this->applySeveralTransforms(DataVector[i]);
        }
    }

    if (file_path != new_path){
        this->eraseFile(new_path);
    }
    
    return DataVector;
}

////////////////////////////////////////////