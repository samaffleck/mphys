add_executable(mphys_tests
    test_main.cpp
    test_field.cpp
    test_mesh.cpp
    test_boundary_condition.cpp
    test_state_vector.cpp
    test_fvm_operators.cpp
    test_model.cpp
    test_solver_options.cpp
    test_sundials_types.cpp
    test_integration_steady.cpp
    test_integration_transient.cpp
    test_spm.cpp
    test_spme.cpp
    test_materials.cpp
    topology_test.cpp
    fvm_mesh_test.cpp
    mesh2d_test.cpp
    mesh_solve_test.cpp
    mesh_coupled_test.cpp
    mesh3d_test.cpp
    mesh_transient_test.cpp
    mesh_precond_test.cpp
    vtk_writer_test.cpp
)

target_link_libraries(mphys_tests
    PRIVATE
    mphys_lib
    GTest::gtest_main
)

# For coverage builds, force every mphys_lib object into the test binary even if
# the current tests don't reference it. Otherwise the linker dead-strips unused
# library code and it never shows up (as uncovered) in the coverage report.
option(MPHYS_COVERAGE "Build for source-based coverage measurement" OFF)
if(MPHYS_COVERAGE)
    target_link_options(mphys_tests PRIVATE
        "-Wl,-force_load,$<TARGET_FILE:mphys_lib>")
endif()

include(GoogleTest)
# Label all mphys tests so they can be run in isolation from SUNDIALS' own
# registered tests:  ctest -L mphys
gtest_discover_tests(mphys_tests PROPERTIES LABELS "mphys")
