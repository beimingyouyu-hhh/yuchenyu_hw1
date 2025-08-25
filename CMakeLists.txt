cmake_minimum_required(VERSION 3.10)

project(main)

# 设置可执行文件输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(algebra src/main.c src/algebra.c)