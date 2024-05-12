cmake_minimum_required(VERSION 3.28)
project(lab15cpp)

set(CMAKE_CXX_STANDARD 17)

add_executable(lab15cpp main.cpp)

include(FetchContent)
FetchContent_Declare(NumCpp
        GIT_REPOSITORY https://github.com/dpilger26/NumCpp
        GIT_TAG Version_2.12.1)
FetchContent_MakeAvailable(NumCpp)

target_link_libraries(${PROJECT_NAME}
        NumCpp::NumCpp
)
