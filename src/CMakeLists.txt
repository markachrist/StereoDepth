cmake_minimum_required(VERSION 2.8)

project(StereoFundus)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

set(UTILS_DIR /Users/mchristopher/Documents/Code/ImageUtils/src)
include_directories(${UTILS_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ -std=c++03")

add_executable(stereo
	main.cpp
	${UTILS_DIR}/RGBUtils.cpp
	DisparityMapper.cpp)

TARGET_LINK_LIBRARIES(stereo ${ITK_LIBRARIES})
