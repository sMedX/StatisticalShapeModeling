cmake_minimum_required(VERSION 3.2)

project(Utils)

set(HEADER
    ssmUtils.h
)

add_custom_target(Utils ALL ${HEADER})
