# include(GNUInstallDirs)

# if(EMSCRIPTEN)
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -sWASM=1 -sASSERTIONS=1 -sFORCE_FILESYSTEM=1")
# endif()

# # Add ImPlot
# include_directories(${CMAKE_SOURCE_DIR}/external/implot)

# # Add GUI include directory
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# # Collect all .cpp files from ImPlot and add them to the build
# file(GLOB IMPlot_SOURCES ${CMAKE_SOURCE_DIR}/external/implot/*.cpp)

# # Collect all .cpp files from src and add them to the build
# file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# # Build the app
# hello_imgui_add_app(StokesFlow StokesFlowMain.cpp ${IMPlot_SOURCES} ${SRC_FILES})

# set_target_properties(StokesFlow PROPERTIES OUTPUT_NAME "AdSim")

# # Link the adsim library
# target_link_libraries(StokesFlow PRIVATE adsim)

# # tinyfiledialogs for native builds
# if(NOT EMSCRIPTEN)
#     add_library(tinyfiledialogs STATIC
#         ${CMAKE_SOURCE_DIR}/external/tinyfiledialogs/tinyfiledialogs.c
#     )
#     target_include_directories(tinyfiledialogs PUBLIC
#         ${CMAKE_SOURCE_DIR}/external/tinyfiledialogs
#     )
#     target_link_libraries(StokesFlow PRIVATE tinyfiledialogs)
# endif()

# # Copy assets next to the built exe (build tree) on Windows
# if(WIN32 AND NOT EMSCRIPTEN AND NOT APPLE)
#     add_custom_command(TARGET StokesFlow POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy_directory
#                 "${CMAKE_CURRENT_SOURCE_DIR}/assets"
#                 "$<TARGET_FILE_DIR:StokesFlow>/assets"
#         COMMENT "Copying assets next to AdSim.exe for local runs")
# endif()

# # Emscripten target flags (if building web)
# if(EMSCRIPTEN)
#     target_compile_options(StokesFlow PRIVATE ${EMSCRIPTEN_FLAGS})
#     target_link_options(StokesFlow PRIVATE ${EMSCRIPTEN_FLAGS})
#     set(ASSET_DIR ${CMAKE_CURRENT_SOURCE_DIR}/assets)
#     set_target_properties(StokesFlow PROPERTIES LINK_FLAGS "--preload-file ${ASSET_DIR}@/assets")
#     target_link_options(
#       StokesFlow PRIVATE
#         "SHELL:-sDISABLE_EXCEPTION_CATCHING=0"
#         "SHELL:-sASSERTIONS=2"
#         "SHELL:-sPTHREAD_POOL_SIZE=4"
#         "SHELL:-sALLOW_MEMORY_GROWTH=1"
#     )
# endif()

# # macOS: copy assets into .app bundle Resources (build tree)
# if(APPLE)
#     set(APP_RES_DIR "${CMAKE_BINARY_DIR}/StokesFlow.app/Contents/Resources")
#     add_custom_command(TARGET StokesFlow POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E make_directory "${APP_RES_DIR}"
#         COMMAND ${CMAKE_COMMAND} -E copy_directory
#                 "${CMAKE_CURRENT_SOURCE_DIR}/assets"
#                 "${APP_RES_DIR}/assets"
#         COMMENT "Copying assets into app bundle Resources"
#     )
# endif()

# # Install the main target
# if(APPLE)
#     install(TARGETS StokesFlow
#         BUNDLE DESTINATION .
#     )
# elseif(EMSCRIPTEN)
#     # usually skip install for web builds
# else()
#     install(TARGETS StokesFlow
#         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#     )
# endif()

# # Install assets
# if(WIN32)
#     install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/assets/
#             DESTINATION ${CMAKE_INSTALL_BINDIR}/assets)
# else()
#     install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/assets
#             DESTINATION ${CMAKE_INSTALL_DATADIR}/AdSim)
# endif()

# # Docs
# install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/README.md ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
#         DESTINATION ${CMAKE_INSTALL_DOCDIR} OPTIONAL)

# # === CPack configuration ===
# include(InstallRequiredSystemLibraries)

# set(CPACK_PACKAGE_NAME "AdSim")
# set(CPACK_PACKAGE_VENDOR "AdSim")
# set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
#     "Adsorption simulation software for modelling pressure and temperature swing adsorption.")
# set(CPACK_PACKAGE_VERSION "2025.5")
# set(CPACK_PACKAGE_INSTALL_DIRECTORY "AdSim")
# set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

# # Generators per platform
# if(WIN32 AND NOT UNIX)
#     set(CPACK_GENERATOR "NSIS")
#     set(CPACK_PACKAGE_EXECUTABLES "AdSim;AdSim")
#     set(CPACK_CREATE_DESKTOP_LINKS "AdSim")
#     set(CPACK_NSIS_MUI_ICON    "${CMAKE_BINARY_DIR}/icon.ico")
#     set(CPACK_NSIS_MUI_UNIICON "${CMAKE_BINARY_DIR}/icon.ico")
#     set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
#     set(CPACK_NSIS_DISPLAY_NAME "AdSim")
# elseif(APPLE)
#     set(CPACK_GENERATOR "DragNDrop")
# elseif(UNIX AND NOT APPLE)
#     set(CPACK_GENERATOR "DEB;RPM;TGZ")
# endif()

# set(CPACK_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION OFF)
# set(CPACK_WARN_ON_ABSOLUTE_INSTALL_DESTINATION ON)

# include(CPack)
