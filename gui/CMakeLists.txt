add_executable(mphys_gui mphys_gui.cpp)

target_link_libraries(mphys_gui PRIVATE mphys_lib implot)

if(APPLE)
    target_link_libraries(mphys_gui PRIVATE "-framework Cocoa" "-framework IOKit" "-framework CoreVideo")
endif()
