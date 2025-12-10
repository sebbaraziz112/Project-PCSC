#include "../includes/DataTypes.hpp"


int computePadding(uint32_t width, uint16_t bitsPerPixel){
    
    int paddingBytes = (4 - ((width * 3) % 4)) % 4;
    return paddingBytes;
}

void initBMPHeader(BMPHeader &header, uint32_t width, uint32_t height, uint16_t bitsPerPixel){
    header.BMP[0] = 'B';
    header.BMP[1] = 'M';
    int padding = computePadding(width, bitsPerPixel);
    header.FileSize = 54 + (width * height * (bitsPerPixel / 8)) + padding * height;
    header.ReservedPartOne = 0;
    header.ReservedPartTwo = 0;
    header.DataOffset = 54;
    header.DIBSize = 40;
    header.Width = width;
    header.Height = height;
    header.Planes = 1;
    header.BitsPerPixel = bitsPerPixel;
    header.Compression = 0;
    header.ImageSize = width * height * (bitsPerPixel / 8) + padding * height;
    header.XPixelsPerMeter = 2835; 
    header.YPixelsPerMeter = 2835; 
    header.ColorsUsed = 0;
    header.ImportantColors = 0;
}
