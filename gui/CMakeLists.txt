add_executable(mphys_gui mphys_gui.cpp)

# tinyfiledialogs — single C file, compiled into the GUI target directly
target_sources(mphys_gui PRIVATE
    ${CMAKE_SOURCE_DIR}/external/tinyfiledialogs/tinyfiledialogs.c
)

set_source_files_properties(
    ${CMAKE_SOURCE_DIR}/external/tinyfiledialogs/tinyfiledialogs.c
    PROPERTIES SKIP_PRECOMPILE_HEADERS ON
)

target_include_directories(mphys_gui PRIVATE
    include
    ${CMAKE_SOURCE_DIR}/external/cereal/include
    ${CMAKE_SOURCE_DIR}/external/tinyfiledialogs
)

target_compile_definitions(mphys_gui PRIVATE
    MPHYS_ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/assets"
)

target_precompile_headers(mphys_gui PRIVATE gui_pch.hpp)

target_link_libraries(mphys_gui PRIVATE mphys_lib implot)

if(APPLE)
    target_link_libraries(mphys_gui PRIVATE "-framework Cocoa" "-framework IOKit" "-framework CoreVideo")
endif()
