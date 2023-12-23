cmake_minimum_required(VERSION 3.12)
project(WATERPXIELS LANGUAGES CXX C)

# Enable C++17
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compilation and optimization flags
if(MSVC)
        SET(CMAKE_CXX_FLAGS_DEBUG "/MP /EHsc /Od /ZI")
        SET(CMAKE_CXX_FLAGS_RELEASE "/MP /EHsc /O2")

        # Set debug dir
        set(CMAKE_VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_ROOT}")
else()
        SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -Wextra")
        SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra")
endif()

file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/include/*.hpp)


file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/third_party/libtim/*)

add_library(libtim ${HEADERS})
target_include_directories(libtim PUBLIC ${CMAKE_SOURCE_DIR}/third_party/libtim)

# For better display in IDEs
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR}/src PREFIX src FILES ${SOURCES})
#SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR}/third_party PREFIX third_party FILES ${LIB})

add_executable(waterpixels ${SOURCES})
target_include_directories(libtim PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_link_libraries(waterpixels PUBLIC libtim)
