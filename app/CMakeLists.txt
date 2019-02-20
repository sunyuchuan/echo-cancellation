cmake_minimum_required(VERSION 2.8)

project(echo_cancellation)
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} ${CMAKE_CXX_FLAGS} -O0 -g2 -ggdb -Wall \
-Wall -Werror -Wno-deprecated-declarations \
-Wno-disabled-macro-expansion -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic \
-Wno-global-constructors -Wno-exit-time-destructors \
-Wno-missing-prototypes -Wno-padded -Wno-old-style-cast")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} ${CMAKE_CXX_FLAGS} -O3 -Wall")
add_definitions(-D __STDC_CONSTANT_MACROS)
add_definitions(-D SAVE_PCM)
set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_C_EXTENSIONS      ON)
set(CMAKE_C_STANDARD        99)
set(CMAKE_CXX_EXTENSIONS    ON)
set(CMAKE_CXX_STANDARD      11)

if(NOT PROJECT_DIR)
    set( PROJECT_DIR ${CMAKE_SOURCE_DIR}/../)
endif()

include_directories(${PROJECT_DIR})

add_executable(${PROJECT_NAME} ${src} 
${PROJECT_DIR}/aec/delay_estimator/delay_estimator_wrapper.c 
${PROJECT_DIR}/aec/delay_estimator/delay_estimator.c 
${PROJECT_DIR}/aec/adp_filter_coeff_factor.c 
${PROJECT_DIR}/aec/aec_control.cc 
${PROJECT_DIR}/aec/aec_core.c 
${PROJECT_DIR}/subband/analy_synth/filterbank.c 
${PROJECT_DIR}/utility/fft/bitrv2.c 
${PROJECT_DIR}/utility/fft/bitrv2conj.c 
${PROJECT_DIR}/utility/fft/cdft_8g.c 
${PROJECT_DIR}/utility/fft/fft_wrapper.c 
${PROJECT_DIR}/utility/fft/ft_ops.c 
${PROJECT_DIR}/utility/fft/rdft_8g_init.c 
${PROJECT_DIR}/utility/fft/rdft_8g_make_table.c 
${PROJECT_DIR}/utility/fft/rdft_8g.c 
${PROJECT_DIR}/utility/ringbuffer/ring_buffer.c 
${PROJECT_DIR}/app/example.cc)
target_link_libraries(${PROJECT_NAME} pthread)
