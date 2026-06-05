cmake_minimum_required(VERSION 3.25)
project(mphys VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_TESTS    "Build tests"    ON)
option(BUILD_EXAMPLES "Build examples" ON)

# Emscripten / WebAssembly configuration.
# C++ exceptions are needed by the core library; memory growth lets the heap
# expand at runtime. We deliberately avoid -pthread so the resulting page can be
# served from any static host (pthreads require SharedArrayBuffer + COOP/COEP).
if(EMSCRIPTEN)
    add_compile_options(-fexceptions -sNO_DISABLE_EXCEPTION_CATCHING=1)
    add_link_options(-fexceptions -sNO_DISABLE_EXCEPTION_CATCHING=1 -sALLOW_MEMORY_GROWTH=1)
endif()

# ============================================================
# Third-party dependencies
# ============================================================

if(EXISTS "${CMAKE_SOURCE_DIR}/external/sundials/CMakeLists.txt")
    # Remove deprecated option names from cache to avoid duplicate-definition warnings
    foreach(_mod ARKODE CVODE CVODES IDA IDAS KINSOL)
        unset(BUILD_${_mod} CACHE)
    endforeach()

    set(SUNDIALS_INDEX_SIZE "64" CACHE STRING "SUNDIALS index type size (32 or 64)")
    set(BUILD_STATIC_LIBS ON CACHE BOOL "")
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "")
    set(EXAMPLES_ENABLE OFF CACHE BOOL "")
    set(EXAMPLES_INSTALL OFF CACHE BOOL "")
    set(SUNDIALS_ENABLE_IDA ON CACHE BOOL "")
    set(SUNDIALS_ENABLE_IDAS ON CACHE BOOL "")
    set(SUNDIALS_ENABLE_CVODE OFF CACHE BOOL "")
    set(SUNDIALS_ENABLE_CVODES OFF CACHE BOOL "")
    set(SUNDIALS_ENABLE_ARKODE OFF CACHE BOOL "")
    set(SUNDIALS_ENABLE_KINSOL ON CACHE BOOL "")

    # Silence warnings in SUNDIALS' own source files
    set(_saved_c_flags "${CMAKE_C_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    add_subdirectory(external/sundials SYSTEM)
    set(CMAKE_C_FLAGS "${_saved_c_flags}")
endif()

# ============================================================
# GUI dependencies (ImGui + ImPlot)
# ============================================================
option(BUILD_GUI "Build ImGui visualisation GUI" ON)

# Suppress warnings in third-party GUI sources across compilers.
if(MSVC)
    set(MPHYS_NO_WARN_FLAG /w)
else()
    set(MPHYS_NO_WARN_FLAG -w)
endif()

if(BUILD_GUI)
    add_library(imgui STATIC
        external/imgui/imgui.cpp
        external/imgui/imgui_draw.cpp
        external/imgui/imgui_tables.cpp
        external/imgui/imgui_widgets.cpp
        external/imgui/backends/imgui_impl_glfw.cpp
        external/imgui/backends/imgui_impl_opengl3.cpp
    )
    target_include_directories(imgui PUBLIC
        external/imgui
        external/imgui/backends
    )

    if(EMSCRIPTEN)
        # Use the contrib.glfw3 port (not the built-in -sUSE_GLFW=3): it has proper
        # HiDPI / canvas-resize support that ImGui_ImplGlfw_InstallEmscriptenCallbacks
        # relies on. WebGL2 gives us a GLES3 context for imgui_impl_opengl3.
        # The flag is needed at compile (for headers) and link (for the JS port).
        target_compile_options(imgui PUBLIC "--use-port=contrib.glfw3")
        target_link_options(imgui PUBLIC "--use-port=contrib.glfw3" "-sMAX_WEBGL_VERSION=2")
    else()
        # Desktop: prefer a system GLFW, otherwise build the vendored submodule.
        find_package(glfw3 QUIET)
        if(NOT glfw3_FOUND)
            if(EXISTS "${CMAKE_SOURCE_DIR}/external/glfw/CMakeLists.txt")
                message(STATUS "System GLFW not found — building vendored external/glfw")
                set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
                set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
                set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
                set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
                # On Linux, build only the X11 backend (works everywhere via XWayland)
                # so we don't drag in the Wayland toolchain (wayland-scanner, protocols).
                # X11 dev headers come from xorg-dev; GLFW_BUILD_X11/WAYLAND are no-ops
                # on non-Linux platforms.
                set(GLFW_BUILD_X11      ON  CACHE BOOL "" FORCE)
                set(GLFW_BUILD_WAYLAND  OFF CACHE BOOL "" FORCE)
                add_subdirectory(external/glfw SYSTEM)
            else()
                message(FATAL_ERROR
                    "GLFW not found and external/glfw submodule is missing.\n"
                    "Run: git submodule update --init external/glfw")
            endif()
        endif()
        find_package(OpenGL REQUIRED)
        target_link_libraries(imgui PUBLIC glfw OpenGL::GL)
    endif()

    target_compile_options(imgui PRIVATE ${MPHYS_NO_WARN_FLAG})

    add_library(implot STATIC
        external/implot/implot.cpp
        external/implot/implot_items.cpp
    )
    target_include_directories(implot PUBLIC external/implot)
    target_link_libraries(implot PUBLIC imgui)
    target_compile_options(implot PRIVATE ${MPHYS_NO_WARN_FLAG})
endif()

# ============================================================
# Core library + executable
# ============================================================
add_subdirectory(src)

# ============================================================
# Examples
# ============================================================
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# ============================================================
# GUI
# ============================================================
if(BUILD_GUI)
    add_subdirectory(gui)
endif()

# ============================================================
# Tests
# ============================================================
if(BUILD_TESTS)
    if(EXISTS "${CMAKE_SOURCE_DIR}/external/googletest/CMakeLists.txt")
        enable_testing()
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        add_subdirectory(external/googletest SYSTEM)
        include(CTest)
        add_subdirectory(tests)
    else()
        message(WARNING "BUILD_TESTS=ON but external/googletest not found — skipping tests")
    endif()
endif()
