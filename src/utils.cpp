#include "../includes/utils.hpp"



bool compareFiles(std::string file1_name, std::string file2_name){
    std::ifstream file1_input(file1_name, std::ios::binary | std::ios::in);
    std::ifstream file2_input(file2_name, std::ios::binary | std::ios::in);

    file1_input.seekg(0);
    file2_input.seekg(0);

    size_t file1_size = file1_input.tellg();
    size_t file2_size = file2_input.tellg();

    if (file1_size != file2_size){
        return false;
    }

    std::vector<char> buffer1(4096), buffer2(4096);

    while (file1_input){
        file1_input.read(buffer1.data(), sizeof(buffer1.data()));
        file2_input.read(buffer2.data(), sizeof(buffer2.data()));

        std::streamsize bytesRead1 = file1_input.gcount(); 
        std::streamsize bytesRead2 = file2_input.gcount(); 

        if (bytesRead1 != bytesRead2) return false; 
        if (!std::equal(buffer1.begin(), buffer1.begin() + bytesRead1, buffer2.begin())) return false;
    }
    return true;
}

void displayPreparation(double upperBound, std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>& matrices){
    const int nChann = matrices.size();
    for (int k = 0; k < nChann; k ++){
        matrices[k] = (matrices[k].array()).log1p().matrix();
        double max = matrices[k].maxCoeff();
        matrices[k] = (matrices[k].array() * (upperBound/max)).matrix();
    }
}
