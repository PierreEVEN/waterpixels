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

file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/src/*.hpp)
file(GLOB_RECURSE LIBTIM ${CMAKE_SOURCE_DIR}/third_party/libtim/*.cpp ${CMAKE_SOURCE_DIR}/third_party/libtim/*.hpp ${CMAKE_SOURCE_DIR}/third_party/libtim/*.h ${CMAKE_SOURCE_DIR}/third_party/libtim/*.hxx)
file(GLOB_RECURSE GLM ${CMAKE_SOURCE_DIR}/third_party/glm/glm/*.cpp ${CMAKE_SOURCE_DIR}/third_party/glm/glm/*.hpp ${CMAKE_SOURCE_DIR}/third_party/glm/glm/*.h ${CMAKE_SOURCE_DIR}/third_party/glm/glm/*.hxx)

# For better display in IDEs
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR}/src PREFIX src FILES ${SOURCES})
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR}/third_party/libtim PREFIX libtim FILES ${LIBTIM})
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR}/third_party/glm PREFIX glm FILES ${GLM})

add_library(glm ${GLM})
target_include_directories(glm PUBLIC ${CMAKE_SOURCE_DIR}/third_party/glm)

add_library(libtim ${LIBTIM})
target_include_directories(libtim PUBLIC ${CMAKE_SOURCE_DIR}/third_party/libtim/)
target_include_directories(libtim PUBLIC ${CMAKE_SOURCE_DIR}/third_party/)

add_executable(waterpixels ${SOURCES})
target_include_directories(waterpixels PRIVATE ${CMAKE_SOURCE_DIR}/src/)
target_link_libraries(waterpixels PRIVATE libtim)
target_link_libraries(waterpixels PRIVATE glm)

set_target_properties(
    waterpixels PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
	