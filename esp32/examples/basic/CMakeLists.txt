# Example project: ARBITER on ESP32 (self-test + optional MPU-9250 over I2C).
# Build:  idf.py set-target esp32 && idf.py build flash monitor
cmake_minimum_required(VERSION 3.16)
set(EXTRA_COMPONENT_DIRS "../..")
include($ENV{IDF_PATH}/tools/cmake/project.cmake)
project(arbiter_basic)
