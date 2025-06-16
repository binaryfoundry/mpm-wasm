cmake_minimum_required(VERSION 3.10)  # newer minimum is recommended

project(mpm-wasm C CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(SRC_DIR "${PROJECT_SOURCE_DIR}")

# Source and header files
set(SOURCES
    src/main.cpp
    src/canvas.cpp
    src/worker.cpp
)

set(HEADERS
    src/main.hpp
    src/canvas.hpp
    src/worker.hpp
)

# Create target
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Organize source groups
source_group("Source" FILES ${SOURCES} ${HEADERS})

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/lib/glm
)

# Add compile options
target_compile_options(${PROJECT_NAME} PRIVATE
    -std=c++11
    -s USE_PTHREADS=1
    -O3
    -pthread
)

# Add link options if building with Emscripten
if (EMSCRIPTEN)
    set_target_properties(
        ${PROJECT_NAME}
        PROPERTIES SUFFIX ".html"
        LINK_FLAGS "--bind -s WASM=1 -s NO_EXIT_RUNTIME=1 -s TOTAL_MEMORY=134217728 -s USE_PTHREADS=1 -s PTHREAD_POOL_SIZE=12 -s ASSERTIONS=2 -s -std=c++11 -O3 --use-preload-plugins --shell-file \"${PROJECT_SOURCE_DIR}/shell.html\"")

    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".html")
endif()
