aux_source_directory(. DIR_SRC)
add_library(chromaxx_clqcd ${DIR_SRC})
target_include_directories(chromaxx_clqcd PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
	$<INSTALL_INTERFACE:include>)
