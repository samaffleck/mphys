add_executable(example_convection_diffusion
    convection_diffusion/main.cpp
)
target_link_libraries(example_convection_diffusion PRIVATE mphys_lib)

add_executable(example_validation_1d_diffusion
    validation_1d_diffusion/main.cpp
)
target_link_libraries(example_validation_1d_diffusion PRIVATE mphys_lib)

add_executable(example_spherical_diffusion
    spherical_diffusion/main.cpp
)
target_link_libraries(example_spherical_diffusion PRIVATE mphys_lib)

add_executable(example_poisson_2d
    poisson_2d/main.cpp
)
target_link_libraries(example_poisson_2d PRIVATE mphys_lib)
