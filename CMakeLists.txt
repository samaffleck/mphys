cmake_minimum_required(VERSION 3.25)
project(mphys VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_TESTS    "Build tests"    ON)
option(BUILD_EXAMPLES "Build examples" ON)

# Emscripten configuration
if(EMSCRIPTEN)
    set(EMSCRIPTEN_FLAGS "-fexceptions" "-sNO_DISABLE_EXCEPTION_CATCHING=1" "-pthread")
    string(JOIN " " EMSCRIPTEN_FLAGS_STR ${EMSCRIPTEN_FLAGS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EMSCRIPTEN_FLAGS_STR}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${EMSCRIPTEN_FLAGS_STR}")
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
