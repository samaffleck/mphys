add_library(mphys_lib STATIC
    mphys.cpp
    # Add further .cpp source files here
)

target_include_directories(mphys_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

if(TARGET SUNDIALS::ida)
    target_link_libraries(mphys_lib PUBLIC SUNDIALS::ida SUNDIALS::idas)
endif()

add_executable(mphys main.cpp)
target_link_libraries(mphys PRIVATE mphys_lib)
