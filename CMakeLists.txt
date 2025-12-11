


cmake_minimum_required (VERSION 3.10)
set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_SHALLOW TRUE    
  GIT_PROGRESS TRUE
)

FetchContent_MakeAvailable(Eigen)


FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()

find_package(Boost REQUIRED COMPONENTS iostreams)


add_library(core 
src/DataReader.cpp
src/DataWriter.cpp
src/FormatFactoring.cpp
src/methods.cpp
src/Transforms.cpp
src/utils.cpp
src/DataTypes.cpp
src/Visual.cpp
)

target_include_directories(core PUBLIC ${CMAKE_SOURCE_DIR}/external/gnuplot-iostream)
target_include_directories(core PUBLIC ${eigen_SOURCE_DIR})


if(Boost_FOUND)
    target_include_directories(core PUBLIC ${Boost_INCLUDE_DIRS})

    target_link_libraries(core
        PUBLIC
            Boost::iostreams
    )
endif()

add_executable(imsonpro main.cpp)
target_link_libraries(imsonpro PRIVATE core)

add_executable(unit_tests tests/data_rw_tests.cpp tests/dft_tests.cpp tests/shift_tests.cpp tests/matrix_vector_tests.cpp)
target_link_libraries(unit_tests PRIVATE core gtest_main)

# tools

add_executable(imsonpro-dft-image tools/dft-image.cpp)
target_link_libraries(imsonpro-dft-image PRIVATE core)

add_executable(imsonpro-convolve-image tools/convolve-image.cpp)
target_link_libraries(imsonpro-convolve-image PRIVATE core)

add_executable(imsonpro-bandfilter-image tools/bandfilter-image.cpp)
target_link_libraries(imsonpro-bandfilter-image PRIVATE core)

add_executable(imsonpro-hist-image tools/histogram-image.cpp)
target_link_libraries(imsonpro-hist-image PRIVATE core)

add_executable(imsonpro-dft-sound tools/dft-sound.cpp)
target_link_libraries(imsonpro-dft-sound PRIVATE core)

add_executable(imsonpro-bandfilter-sound tools/bandfilter-sound.cpp)
target_link_libraries(imsonpro-bandfilter-sound PRIVATE core)

add_executable(imsonpro-hist-sound tools/histogram-sound.cpp)
target_link_libraries(imsonpro-hist-sound PRIVATE core)

include(GoogleTest)
gtest_discover_tests(unit_tests)

