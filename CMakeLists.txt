cmake_minimum_required (VERSION 2.6)
project (score_meta)

#FIND_PACKAGE( Boost 1.46 COMPONENTS program_options system log log_setup thread REQUIRED )
FIND_PACKAGE( Boost 1.46 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

list( APPEND CMAKE_CXX_FLAGS "-std=gnu++11 -O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
#list( APPEND CMAKE_CXX_FLAGS "-std=gnu++11 -g -O3 -pg -pedantic-errors ${CMAKE_CXX_FLAGS}")

add_executable(score_meta score_meta.cpp)

add_library(Seed Seed.cpp)
add_library(OneSequence OneSequence.cpp)
add_library(AllSequences AllSequences.cpp)
target_link_libraries(score_meta AllSequences OneSequence Seed ${Boost_LIBRARIES})
