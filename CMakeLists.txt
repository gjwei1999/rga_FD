cmake_minimum_required(VERSION 3.5)
project(clas12_analysis)

set(CMAKE_CXX_STANDARD 14)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
include(${ROOT_USE_FILE})

include_directories(src/include)

file(GLOB clas12lib_SRC
    "src/*.cpp"
)
add_library(clas12lib STATIC ${clas12lib_SRC})


file(GLOB clas12_EXE
    "exe/*.cpp"
)

foreach(exe_file IN LISTS clas12_EXE)
    get_filename_component(exe_name ${exe_file} NAME_WE)
    add_executable(${exe_name} ${exe_file})
    target_link_libraries(${exe_name} clas12lib ${ROOT_LIBRARIES})
endforeach()
