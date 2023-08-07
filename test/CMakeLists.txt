cmake_minimum_required(VERSION 3.10)
project(EmbeddedQuaternion_Demo)

set(CMAKE_C_STANDARD 11)

# 指定生成目标
add_executable(Quaternion_test Quaternion_test.c ../src/Quaternion.c)
target_link_libraries(Quaternion_test PRIVATE m)
