cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME mpm-wasm)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(SRC_DIR "${PROJECT_SOURCE_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_STANDARD 11)

if (EMSCRIPTEN)
    add_definitions(-O3)
endif ()

project(${PROJECT_NAME} C CXX)

set(SOURCES
    src/main.cpp
    src/canvas.cpp
    src/worker.cpp)

set(HEADERS
    src/main.hpp
    src/canvas.hpp
    src/worker.hpp)

SOURCE_GROUP("Source" FILES ${SOURCES})
SOURCE_GROUP("Source" FILES ${HEADERS})

include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/lib/glm)

add_executable(
    ${PROJECT_NAME}
    ${SOURCES}
    ${HEADERS})

target_link_libraries(
    ${PROJECT_NAME})

if (EMSCRIPTEN)
    set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES SUFFIX ".html"
        LINK_FLAGS "--bind -s WASM=1 -s NO_EXIT_RUNTIME=1 -s TOTAL_MEMORY=134217728 -s USE_PTHREADS=1 -s PTHREAD_POOL_SIZE=12 -s ASSERTIONS=2 -s --emrun -std=c++11 -O3 --use-preload-plugins --shell-file \"../../shell.html\"")
endif ()
