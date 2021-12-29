# listing the TENDON-DRIVEN library source files into a variable
set(TENDONDRIVEN_LIB_SOURCES
    src/TendonDriven.cpp
    src/ODESystem.cpp
    src/PDESystem.cpp
    src/timeManager.cpp
)

# creating a TENDON-DRIVEN static library
add_library(TENDONDRIVENlib STATIC ${TENDONDRIVEN_LIB_SOURCES})
target_include_directories(TENDONDRIVENlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)