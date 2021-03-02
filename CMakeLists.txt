# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 3.18)

# This project is designed to be built outside the Insight source tree.
project(IsoTest)

# Find ITK.
#set(ITK_DIR, C:\\InsightToolkit-5.1.2\\bin)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_library(IsoLib STATIC 
			R3Graph.h R3Graph.cpp
			voxelset.h voxelset.cpp
			ddd.cpp
			roi.h roi.cpp
			pi.h
			Triangulation.h Triangulation.cpp
			matrix.h
			r2geom.h r2geom.cpp
			skala.h skala.cpp
			BinHeap.h BinHeap.cpp)

add_executable(IsoTest startUp.cxx)

target_link_libraries(IsoTest ${ITK_LIBRARIES}  IsoLib)
