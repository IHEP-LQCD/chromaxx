# basic setup for cmake
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

execute_process(COMMAND chroma-config --cxx
OUTPUT_VARIABLE CHROMA_CXX OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND chroma-config --ldflags
OUTPUT_VARIABLE CHROMA_LDFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND chroma-config --cxxflags
OUTPUT_VARIABLE CHROMA_CXXFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND chroma-config --libs
OUTPUT_VARIABLE CHROMA_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

message("CHROMA_CXX ${CHROMA_CXX}")
message("CHROMA_LDFLAGS ${CHROMA_LDFLAGS}")
message("CHROMA_CXXFLAGS ${CHROMA_CXXFLAGS}")
message("CHROMA_LIBS ${CHROMA_LIBS}")

PROJECT(chroma_addons)				            #设置工程名
MESSAGE(STATUS "Project chroma_addons")	 	    #打印消息

set(CMAKE_CXX_FLAGS ${CHROMA_CXXFLAGS})

include_directories("./include")

add_subdirectory(src)
add_executable(chroma_addons chroma.cc register_new.cc)
target_link_libraries(chroma_addons chroma_clqcd ${CHROMA_LDFLAGS} ${CHROMA_LIBS})
add_executable(purgauge_addons purgaug.cc register_new.cc)
target_link_libraries(purgauge_addons chroma_clqcd ${CHROMA_LDFLAGS} ${CHROMA_LIBS})
install(TARGETS chroma_addons)
