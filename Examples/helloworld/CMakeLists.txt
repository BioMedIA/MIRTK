cmake_minimum_required(VERSION 3.9)
project(flip)

find_package(VTK)

find_package(MIRTK REQUIRED
	PointSet
	Common 
	Numerics 
	IO 
	Image 
	Transformation 
	Registration
	)

include(MIRTKTargets)

add_executable(flip flip.cc)
target_include_directories(flip PUBLIC "${MIRTK_INCLUDE_DIRS}")
target_link_libraries(flip mirtk::LibCommon mirtk::LibPointSet)
install(TARGETS flip DESTINATION bin)
