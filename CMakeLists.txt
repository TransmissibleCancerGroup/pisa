project(pisa)
cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 11)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/deps/SeqLib/src/svaba/SeqLib/src/libseqlib.a")
    message(${CMAKE_SOURCE_DIR}/build_deps.sh " " $ENV{CPPFLAGS})
    execute_process(COMMAND ${CMAKE_SOURCE_DIR}/build_deps.sh " " $ENV{CPPFLAGS})
endif()

find_package(Threads REQUIRED)

include_directories(
        deps/SeqLib
        deps/SeqLib/htslib
        deps/CTPL
        src      
)

set(DEPLIBS
        ${CMAKE_SOURCE_DIR}/deps/SeqLib/src/libseqlib.a
        ${CMAKE_SOURCE_DIR}/deps/SeqLib/bwa/libbwa.a
        ${CMAKE_SOURCE_DIR}/deps/SeqLib/htslib/libhts.a
        ${CMAKE_SOURCE_DIR}/deps/SeqLib/fermi-lite/libfml.a)

set(SOURCE_FILES src/main.cpp src/readCorrection.cpp src/readFilter.cpp src/threadedFermiAssembler.cpp)
add_executable(pisa ${SOURCE_FILES})
target_link_libraries(pisa z lzma bz2 Threads::Threads ${DEPLIBS})