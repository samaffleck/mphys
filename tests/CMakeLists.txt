add_executable(mphys_tests
    test_main.cpp
    topology_test.cpp
    fvm_mesh_test.cpp
)

target_link_libraries(mphys_tests
    PRIVATE
    mphys_lib
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(mphys_tests)
