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
  gnuplot_iostream
  GIT_REPOSITORY https://github.com/dstahlke/gnuplot-iostream.git
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(gnuplot_iostream)

# -------------------- Boost --------------------
# On Linux : STATIC ou SHARED selon ton installation
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)


FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()



add_library(core 
src/DataReader.cpp
src/DataWriter.cpp
src/FormatFactoring.cpp
src/methods.cpp
src/Transforms.cpp
src/utils.cpp
src/Visual.cpp
)

target_include_directories(core PUBLIC ${eigen_SOURCE_DIR} ${gnuplot_iostream_SOURCE_DIR}
${Boost_INCLUDE_DIRS}
)

# Linker Boost sur la librairie core
target_link_libraries(core PRIVATE 
    ${Boost_IOSTREAMS_LIBRARY} 
    ${Boost_SYSTEM_LIBRARY} 
    ${Boost_FILESYSTEM_LIBRARY}
)


add_executable(imsonpro main.cpp)
target_link_libraries(imsonpro PRIVATE core)

add_executable(unit_tests tests/data_rw_tests.cpp tests/dft_tests.cpp tests/shift_tests.cpp)
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

find_package(Doxygen REQUIRED)

set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/docs/Doxyfile.in)
set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

set(DOC_HTML_DIR ${CMAKE_BINARY_DIR}/html)

add_custom_command(
    OUTPUT ${DOC_HTML_DIR}/index.html
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Génération de la documentation Doxygen..."
    VERBATIM
)

add_custom_target(doc ALL
    DEPENDS ${DOC_HTML_DIR}/index.html
    COMMENT "Documentation Doxygen prête → ${DOC_HTML_DIR}/index.html"
)

# Ouverture automatique pour WSL (Linux)
add_custom_command(TARGET doc POST_BUILD
    COMMAND wslview ${DOC_HTML_DIR}/index.html
    COMMENT "Ouverture automatique de la documentation dans Windows (via WSL)..."
)
