add_executable(mphys_gui mphys_gui.cpp)

target_include_directories(mphys_gui PRIVATE
    include
    ${CMAKE_SOURCE_DIR}/external/cereal/include
)

target_precompile_headers(mphys_gui PRIVATE gui_pch.hpp)

target_link_libraries(mphys_gui PRIVATE mphys_lib implot)

if(EMSCRIPTEN)
    # ---- WebAssembly build ----------------------------------------------------
    # Emit a runnable HTML shell, bundle the assets into the virtual filesystem at
    # /assets, and let the canvas resize with the browser window.
    set_target_properties(mphys_gui PROPERTIES SUFFIX ".html")
    target_compile_definitions(mphys_gui PRIVATE MPHYS_ASSETS_DIR="/assets")
    target_link_options(mphys_gui PRIVATE
        "--preload-file=${CMAKE_CURRENT_SOURCE_DIR}/assets@/assets"
        "--shell-file=${CMAKE_CURRENT_SOURCE_DIR}/shell.html"
        "-sALLOW_MEMORY_GROWTH=1"
    )
    # Rebuild the page if the shell template changes.
    set_property(TARGET mphys_gui APPEND PROPERTY LINK_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/shell.html")
else()
    # ---- Desktop build (macOS / Windows / Linux) ------------------------------
    # tinyfiledialogs provides native Open/Save dialogs; it is a single C file
    # compiled straight into the GUI target (no web equivalent).
    target_sources(mphys_gui PRIVATE
        ${CMAKE_SOURCE_DIR}/external/tinyfiledialogs/tinyfiledialogs.c
    )
    set_source_files_properties(
        ${CMAKE_SOURCE_DIR}/external/tinyfiledialogs/tinyfiledialogs.c
        PROPERTIES SKIP_PRECOMPILE_HEADERS ON
    )
    target_include_directories(mphys_gui PRIVATE
        ${CMAKE_SOURCE_DIR}/external/tinyfiledialogs
    )
    target_compile_definitions(mphys_gui PRIVATE
        MPHYS_ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/assets"
    )

    if(APPLE)
        target_link_libraries(mphys_gui PRIVATE "-framework Cocoa" "-framework IOKit" "-framework CoreVideo")
    endif()
endif()
