add_library(mphys_lib STATIC
    solver_options.cpp
    mesh.cpp
    topology.cpp
    state_vector.cpp
    model.cpp
    fvm_operators.cpp
    fvm_mesh.cpp
    mesh_model.cpp
    mesh_steady_solver.cpp
    mesh_transient_solver.cpp
    vtk_writer.cpp
    transient_solver.cpp
    steady_solver.cpp
)

target_include_directories(mphys_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

# Namespaced alias so consumers can link `mphys::mphys` whether mphys is pulled
# in via FetchContent / add_subdirectory or (eventually) find_package.
add_library(mphys::mphys ALIAS mphys_lib)

target_precompile_headers(mphys_lib PRIVATE pch.hpp)

if(TARGET SUNDIALS::ida)
    target_link_libraries(mphys_lib PUBLIC SUNDIALS::ida SUNDIALS::idas)
endif()

if(TARGET SUNDIALS::kinsol)
    target_link_libraries(mphys_lib PUBLIC SUNDIALS::kinsol)
endif()

add_executable(mphys main.cpp)
target_link_libraries(mphys PRIVATE mphys_lib)
