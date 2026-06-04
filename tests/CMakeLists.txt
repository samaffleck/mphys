add_executable(mphys_tests
    test_main.cpp
    test_spm.cpp
)

target_link_libraries(mphys_tests
    PRIVATE
    mphys_lib
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(mphys_tests)
