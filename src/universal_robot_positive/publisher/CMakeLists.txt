cmake_minimum_required(VERSION 2.8.3)
project(publisher)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  kdl_parser
  message_generation

)



catkin_package(
   CATKIN_DEPENDS message_runtime
   )



# check c++11 / c++0x
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "-std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "-std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler. Suggested solution: update the pkg build-essential ")
endif()




file(GLOB folder_source src/*.c src/*.cpp src/*.h )
## source_group("Source Files" FILES ${folder_source})

SET(qpOASES_INCLUDE_DIRS /home/w/qpOASES-3.2.0/include)
SET(qpOASES_LIBRARIES /home/w/qpOASES-3.2.0/bin)

SET(CGAL_INCLUDE_DIRS /usr/local/include)
SET(CGAL_LIBRARIES /usr/local/lib)


include_directories(include ${catkin_INCLUDE_DIRS}
                            ${qpOASES_INCLUDE_DIRS}
                            ${CGAL_INCLUDE_DIRS})

LINK_DIRECTORIES(${qpOASES_LIBRARIES})
LINK_DIRECTORIES(${CGAL_LIBRARIES})

add_executable(talker ${folder_source})




target_link_libraries(talker ${catkin_LIBRARIES} qpOASES CGAL CGAL_Core CGAL_ImageIO gmp mpfr)

add_dependencies(talker ${catkin_EXPORTED_TARGETS})


