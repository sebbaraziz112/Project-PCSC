#ifndef DATATYPES
#define DATATYPES

#include <cstdint>
#include <cmath>

/// @file DataTypes.hpp
/// @brief Defines low-level binary structures for WAV, BMP and simple point formats.


// ============================================================================
// WAV HEADER
// ============================================================================

/**
 * @brief Binary representation of a standard WAV file header.
 *
 * All fields are packed to match the exact byte layout of PCM WAV files.
 * This structure is used for parsing WAV headers and retrieving metadata
 * such as sample rate, channels, and bit depth.
 *
 * @note `#pragma pack(push, 1)` ensures no padding is added by the compiler.
 */
#pragma pack(push, 1)
struct WAVHeader {

    char Riff[4];            ///< Contains the characters "RIFF".
    uint32_t FileSize;       ///< Total file size minus 8 bytes.
    char FileType[4];        ///< Contains the characters "WAVE".
    char Fmt[4];             ///< Contains the characters "fmt ".
    uint32_t FmtChunkSize;   ///< Size of the fmt chunk (typically 16 for PCM).
    uint16_t PCM;            ///< Audio format code (1 = PCM).
    uint16_t NumChannels;    ///< Number of audio channels.
    uint32_t SampleRate;     ///< Samples per second.
    uint32_t ByteRate;       ///< Bytes per second (SampleRate * NumChannels * BitsPerSample/8).
    uint16_t BlockAlign;     ///< Bytes per sample frame.
    uint16_t BitsPerSample;  ///< Bits per audio sample (8, 16, 24, 32...).
    char Data[4];            ///< Contains the characters "data".
    uint32_t DataSize;       ///< Size of the audio data that follows.
};
#pragma pack(pop)


// ============================================================================
// BMP HEADER
// ============================================================================

/**
 * @brief Binary representation of a BMP file header (Bitmap File Header + DIB Header).
 *
 * Contains metadata describing the BMP file, including image dimensions,
 * pixel depth, compression, and size. This structure is used for raw BMP parsing.
 *
 * @warning Only supports the standard BITMAPINFOHEADER format (40 bytes).
 */
#pragma pack(push, 1)
struct BMPHeader {

    // --- BITMAP FILE HEADER ---
    char BMP[2];                            ///< Magic number: must be 'B' 'M'.
    uint32_t FileSize;                      ///< Total file size in bytes.
    uint16_t ReservedPartOne;               ///< Reserved; must be zero.
    uint16_t ReservedPartTwo;               ///< Reserved; must be zero.
    uint32_t DataOffset;                    ///< Offset from start of file to pixel array.

    // --- BITMAP INFO HEADER (DIB HEADER) ---
    uint32_t DIBSize;                       ///< Size of DIB header (40 bytes for BITMAPINFOHEADER).
    uint32_t Width;                         ///< Bitmap width in pixels.
    uint32_t Height;                        ///< Bitmap height in pixels.
    uint16_t Planes;                        ///< Number of color planes (must be 1).
    uint16_t BitsPerPixel;                  ///< Bits per pixel (e.g., 24 or 32).
    uint32_t Compression;                   ///< Compression type (0 = no compression).
    uint32_t ImageSize;                     ///< Raw bitmap data size (may be 0 for uncompressed BMP).
    uint32_t XPixelsPerMeter;               ///< Horizontal resolution (pixels per meter).
    uint32_t YPixelsPerMeter;               ///< Vertical resolution (pixels per meter).
    uint32_t ColorsUsed;                    ///< Number of colors in the palette (0 = default).
    uint32_t ImportantColors;               ///< Number of important colors (0 = all).
};
#pragma pack(pop)


// ============================================================================
// POINT STRUCTURES
// ============================================================================

/**
 * @brief Represents a 2D point with 8-bit unsigned integer coordinates.
 */
#pragma pack(push, 1)
struct Point2Duint8_t {
    uint8_t x;   ///< X coordinate.
    uint8_t y;   ///< Y coordinate.
};
#pragma pack(pop)

/**
 * @brief Represents a 2D point with 16-bit unsigned integer coordinates.
 */
#pragma pack(push, 1)
struct Point2Duint16_t {
    uint16_t x;  ///< X coordinate.
    uint16_t y;  ///< Y coordinate.
};
#pragma pack(pop)

/**
 * @brief Represents a 3D point with 8-bit unsigned integer coordinates.
 */
#pragma pack(push, 1)
struct Point3Duint8_t {
    uint8_t x;   ///< X coordinate.
    uint8_t y;   ///< Y coordinate.
    uint8_t z;   ///< Z coordinate.
};
#pragma pack(pop)

/**
 * @brief Represents a 3D point with 16-bit unsigned integer coordinates.
 */
#pragma pack(push, 1)
struct Point3Duint16_t {
    uint16_t x;  ///< X coordinate.
    uint16_t y;  ///< Y coordinate.
    uint16_t z;  ///< Z coordinate.
};
#pragma pack(pop)

#endif
