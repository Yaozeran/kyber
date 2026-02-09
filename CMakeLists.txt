cmake_minimum_required(VERSION 3.20)

project(kyber)

add_library(kyber INTERFACE)

target_include_directories(kyber INTERFACE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  
)

target_compile_features(kyber INTERFACE
    cxx_std_20
)

add_executable(main test/main.cc)
target_link_libraries(main PRIVATE
  kyber
)

add_executable(demo test/demo.cc)
target_link_libraries(demo PRIVATE
  kyber
)