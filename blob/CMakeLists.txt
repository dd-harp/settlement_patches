cmake_minimum_required(VERSION 3.13)  # Maybe 3.14 required?
project(blob)

set(CMAKE_CXX_STANDARD 14)

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

find_library(
        LIBTIFF
        NAMES tiff
        PATHS /usr/lib/x86_64-linux-gnu
        )

add_executable(blob main.cpp read_tiff.h connected_settlements.h pixel.h)
add_to_cached_list(CGAL_EXECUTABLE_TARGETS blob )

# Link the executable to CGAL and third-party libraries
target_link_libraries(blob   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} tiff )
