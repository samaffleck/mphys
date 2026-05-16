add_executable(mphys_gui mphys_gui.cpp)

target_include_directories(mphys_gui PRIVATE include)

target_compile_definitions(mphys_gui PRIVATE
    MPHYS_ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/assets"
)

target_precompile_headers(mphys_gui PRIVATE gui_pch.hpp)

target_link_libraries(mphys_gui PRIVATE mphys_lib implot)

if(APPLE)
    target_link_libraries(mphys_gui PRIVATE "-framework Cocoa" "-framework IOKit" "-framework CoreVideo")
endif()
