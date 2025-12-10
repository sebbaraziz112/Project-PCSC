#ifndef DATATYPES
#define DATATYPES


#include <cstdint>
#include <cmath>

// WAV Header Structure

#pragma pack(push, 1)
struct WAVHeader {

    char Riff[4];            // 4 bytes
    uint32_t FileSize;       // 4 bytes
    char FileType[4];        // 4 bytes
    char Fmt[4];             // 4 bytes
    uint32_t FmtChunkSize;   // 4 bytes
    uint16_t PCM;            // 2 bytes
    uint16_t NumChannels;    // 2 bytes
    uint32_t SampleRate;     // 4 bytes
    uint32_t ByteRate;       // 4 bytes
    uint16_t BlockAlign;     // 2 bytes
    uint16_t BitsPerSample;  // 2 bytes
    char Data[4];            // 4 bytes
    uint32_t DataSize;       // 4 bytes

};
#pragma pack(pop)

// BMP Header Structure

#pragma pack(push, 1)
struct BMPHeader {

    // First Part: BITMAP FILE HEADER
    char BMP[2];                            // 2 bytes
    uint32_t FileSize;                      // 4 bytes 
    uint16_t ReservedPartOne;               // 2 bytes
    uint16_t ReservedPartTwo;               // 2 bytes
    uint32_t DataOffset;                    // 4 bytes

    // Second Part: BITMAP INFO HEADER (DIB HEADER)
    uint32_t DIBSize;                       // 4 bytes
    uint32_t Width;                         // 4 bytes
    uint32_t Height;                        // 4 bytes
    uint16_t Planes;                        // 2 bytes
    uint16_t BitsPerPixel;                  // 2 bytes
    uint32_t Compression;                   // 4 bytes
    uint32_t ImageSize;                     // 4 bytes
    uint32_t XPixelsPerMeter;               // 4 bytes
    uint32_t YPixelsPerMeter;               // 4 bytes
    uint32_t ColorsUsed;                    // 4 bytes
    uint32_t ImportantColors;               // 4 bytes

};
#pragma pack(pop)

// Some Points

#pragma pack(push, 1)
struct Point2Duint8_t {

    uint8_t x;
    uint8_t y;

};
#pragma pack(pop)

#pragma pack(push, 1)
struct Point2Duint16_t {

    uint16_t x;
    uint16_t y;

};
#pragma pack(pop)

#pragma pack(push, 1)
struct Point3Duint8_t {

    uint8_t x;
    uint8_t y;
    uint8_t z;


};
#pragma pack(pop)

#pragma pack(push, 1)
struct Point3Duint16_t {

    uint16_t x;
    uint16_t y;
    uint16_t z;


};
#pragma pack(pop)

#endif