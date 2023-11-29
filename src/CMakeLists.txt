aux_source_directory(. DIR_SRC)
add_library(chromaxx_clqcd ${DIR_SRC})
set_target_properties(chromaxx_clqcd PROPERTIES CXX_STANDARD "20")
target_link_libraries(chromaxx_clqcd PUBLIC Chroma::chromalib)
target_include_directories(chromaxx_clqcd PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
	$<INSTALL_INTERFACE:include>)
