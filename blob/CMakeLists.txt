cmake_minimum_required(VERSION 3.13)  # Maybe 3.14 required?
project(blob)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
message(STATUS "Cmake module path ${CMAKE_MODULE_PATH}")
include(GoogleTest)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif()



# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

    message(STATUS "This project requires the Boost library, and will not be compiled.")

    return()

endif()

# Geographic lib
find_package (GeographicLib 1.34 REQUIRED)
include_directories (${GeographicLib_INCLUDE_DIRS})
add_definitions (${GeographicLib_DEFINITIONS})


find_library(
        LIBTIFF
        NAMES tiff
        PATHS /usr/lib/x86_64-linux-gnu
        )

add_executable(blob main.cpp read_tiff.h connected_settlements.h pixel.h simple_patches.h
        shapefile.h gdal_raster.h gdal_raster.cpp gdal_vector.h simple_patches_test.cpp
        gdal_vector.cpp projection.h projection.cpp admin_patch.h admin_patch.cpp
        on_demand_raster.h on_demand_raster.cpp sparse_settlements.cpp sparse_settlements.h
        split_patches.h split_patches.cpp metis_io.cpp metis_io.h
        )
add_to_cached_list(CGAL_EXECUTABLE_TARGETS blob )

# Link the executable to CGAL and third-party libraries
target_link_libraries(
        blob
        ${CGAL_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
        ${GeographicLib_LIBRARIES}
        geotiff tiff gdal
        boost_system boost_program_options boost_filesystem
        gtest_main
)

# Make a second executable that runs the tests.
add_executable(blob_test test_main.cpp cgal_tutorial.cpp simple_patches_test.cpp
        simple_patches.h admin_patch_test.cpp admin_patch.cpp admin_patch.h
        projection.h projection.cpp gdal_raster.cpp projection_test.cpp
        on_demand_raster.h on_demand_raster.cpp sparse_settlements.cpp
        sparse_settlements.h split_patches.h split_patches.cpp
        sparse_settlements_test.cpp gdal_raster_test.cpp metis_io.cpp)
target_link_libraries(
        blob_test
        ${CGAL_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
        geotiff tiff gdal
        boost_system boost_program_options boost_filesystem
        gtest_main
)

gtest_discover_tests(blob_test)
