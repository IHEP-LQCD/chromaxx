# basic setup for cmake
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(chromaxx VERSION 1.0.0 LANGUAGES CXX C)

# need Chroma package
find_package(Chroma REQUIRED)

MESSAGE(STATUS "Project chromaxx")

add_subdirectory(src)
add_executable(chromaxx chroma.cc register_new.cc)
#target_link_libraries(chromaxx chromaxx_clqcd ${CHROMA_LDFLAGS} ${CHROMA_LIBS})
set_target_properties(chromaxx PROPERTIES CXX_STANDARD "20")
target_link_libraries(chromaxx chromaxx_clqcd)
add_executable(purgaugexx purgaug.cc register_new.cc)
set_target_properties(purgaugexx PROPERTIES CXX_STANDARD "20")
target_link_libraries(purgaugexx chromaxx_clqcd)
add_executable(hmcxx hmc.cc register_new.cc)
set_target_properties(hmcxx PROPERTIES CXX_STANDARD "20")
target_link_libraries(hmcxx chromaxx_clqcd)
install(TARGETS chromaxx purgaugexx hmcxx chromaxx_clqcd
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib)
install(DIRECTORY include/ DESTINATION include
		FILES_MATCHING PATTERN "*.h")
