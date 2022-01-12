
find_package(Boost COMPONENTS filesystem system iostreams REQUIRED)

# Check for the Boost libray, if found print message, include dirs and link libraries.
if (Boost_FOUND)
    message(STATUS "Boost has been found!")
    include_directories(${Boost_INCLUDE_DIRS})
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/gnuplot_lib")
    link_libraries(${Boost_LIBRARIES})

    if (WIN32)
        # disable autolinking in boost
        add_definitions( -DBOOST_ALL_NO_LIB )  
        # asks Boost to report information about automatic linking requests on Visual Studio
        add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})
        # force all boost libraries to dynamic link
        add_definitions( -DBOOST_ALL_DYN_LINK )
    endif()
elseif(NOT Boost_FOUND)
    error(STATUS "Boost has not been found!")
endif()


# Target.
add_library(gnuplot_iostream INTERFACE)
target_include_directories(gnuplot_iostream INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    $<INSTALL_INTERFACE:>)
target_link_libraries(gnuplot_iostream INTERFACE
    Boost::iostreams
    Boost::system
    Boost::filesystem
)