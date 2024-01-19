cmake_minimum_required(VERSION 3.12)
project(Waterpixels LANGUAGES CXX C)

# Enable C++17
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compilation and optimization flags
if(MSVC)
        SET(CMAKE_CXX_FLAGS_DEBUG "/MP /EHsc /Od /ZI")
        SET(CMAKE_CXX_FLAGS_RELEASE "/MP /EHsc /O2")
else()
        SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -Wextra")
        SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra")
endif()

file(GLOB_RECURSE MAIN ${CMAKE_SOURCE_DIR}/src/main.cpp)
file(GLOB_RECURSE WATERPIXELS ${CMAKE_SOURCE_DIR}/src/waterpixels/*.cpp)
file(GLOB_RECURSE WATERPIXELS_PUBLIC ${CMAKE_SOURCE_DIR}/include/waterpixels/*.hpp)
file(GLOB_RECURSE LIBTIM ${CMAKE_SOURCE_DIR}/third_party/libtim/*.cpp ${CMAKE_SOURCE_DIR}/third_party/libtim/*.hpp ${CMAKE_SOURCE_DIR}/third_party/libtim/*.h ${CMAKE_SOURCE_DIR}/third_party/libtim/*.hxx)
file(GLOB_RECURSE GLM ${CMAKE_SOURCE_DIR}/third_party/glm/glm/*.cpp ${CMAKE_SOURCE_DIR}/third_party/glm/glm/*.hpp ${CMAKE_SOURCE_DIR}/third_party/glm/glm/*.h ${CMAKE_SOURCE_DIR}/third_party/glm/glm/*.hxx)

# For better display in IDEs
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR}/src PREFIX src FILES ${WATERPIXELS})
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR}/include PREFIX include FILES ${WATERPIXELS_PUBLIC})
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR}/third_party/libtim PREFIX libtim FILES ${LIBTIM})
SOURCE_GROUP(TREE ${CMAKE_SOURCE_DIR}/third_party/glm PREFIX glm FILES ${GLM})

add_library(glm ${GLM})
target_include_directories(glm PUBLIC ${CMAKE_SOURCE_DIR}/third_party/glm)

add_library(libtim ${LIBTIM})
target_include_directories(libtim PUBLIC ${CMAKE_SOURCE_DIR}/third_party/libtim/)
target_include_directories(libtim PUBLIC ${CMAKE_SOURCE_DIR}/third_party/)

find_package(OpenMP)
if (OPENMP_FOUND)
	if (MSVC)
        SET(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG} " /openmp:llvm")
        SET(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} " /openmp:llvm")
	else()
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
	message("-- enabled OpenMP")
endif()

add_library(waterpixels ${WATERPIXELS} ${WATERPIXELS_PUBLIC})
target_include_directories(waterpixels PUBLIC ${CMAKE_SOURCE_DIR}/include/)
target_link_libraries(waterpixels PUBLIC libtim)
target_link_libraries(waterpixels PUBLIC glm)

add_executable(main ${MAIN})
target_link_libraries(main PRIVATE waterpixels)

set_target_properties(main PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
	