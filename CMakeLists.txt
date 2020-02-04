cmake_minimum_required(VERSION 3.12)
project(SVD_OMP)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# Compiler Abhängigkeiten
#if(CMAKE_SIZEOF_VOID_P EQUAL 4)
#    message(FATAL_ERROR "32 bit architecture is not supported")
#endif()

#set(CID ${CMAKE_CXX_COMPILER_ID})
#if(WIN32 AND NOT "${CID}" STREQUAL "GNU")
#    if("${CID}" STREQUAL "Intel")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:native /fp:fast=2")
#    elseif("${CID}" STREQUAL "MSVC")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2 /fp:fast")
#    else()
#        message(FATAL_ERROR "Provide your own compile flags.")
#    endif()

#    if("${CID}" STREQUAL "Intel" OR "${CID}" STREQUAL "MSVC")
#        foreach(flag_var CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_RELEASE
#                CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
#            if(${flag_var} MATCHES "/MD")
#                string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
#            endif()
#        endforeach()
#    endif()
#else()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -march=native -fsanitize=address")
#    if("${CID}" STREQUAL "GNU")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -fsanitize=address")
#    endif()
#endif()

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fsanitize=address")

# RUN
add_executable(SVD_OMP_RUN Main/svd.cpp Main/svd.h Main/Main.cpp)
target_link_libraries(SVD_OMP_RUN PUBLIC OpenMP::OpenMP_CXX)

# Test
add_executable(SVD_OMP_TEST Test/test.cpp Main/svd.cpp Main/svd.h Test/catch_main.cpp)
add_subdirectory(catch)
target_link_libraries(SVD_OMP_TEST PRIVATE catch)
