add_library(mphys_lib STATIC
    solver_options.cpp
    mesh.cpp
    state_vector.cpp
    model.cpp
    fvm_operators.cpp
    transient_solver.cpp
    steady_solver.cpp
)

target_include_directories(mphys_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_precompile_headers(mphys_lib PRIVATE pch.hpp)

if(TARGET SUNDIALS::ida)
    target_link_libraries(mphys_lib PUBLIC SUNDIALS::ida SUNDIALS::idas)
endif()

if(TARGET SUNDIALS::kinsol)
    target_link_libraries(mphys_lib PUBLIC SUNDIALS::kinsol)
endif()

add_executable(mphys main.cpp)
target_link_libraries(mphys PRIVATE mphys_lib)
