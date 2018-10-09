cmake_minimum_required(VERSION 2.8)

project(echo_cancellation)
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} ${CMAKE_CXX_FLAGS} -O0 -g2 -ggdb -Wall \
-Wno-deprecated-declarations \
-Wno-disabled-macro-expansion -Wno-float-equal -Wno-c++98-compat -Wno-c++98-compat-pedantic \
-Wno-global-constructors -Wno-exit-time-destructors \
-Wno-missing-prototypes -Wno-padded -Wno-old-style-cast")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} ${CMAKE_CXX_FLAGS} -O3 -Wall")
add_definitions(-D __STDC_CONSTANT_MACROS)
set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_C_EXTENSIONS      ON)
set(CMAKE_C_STANDARD        99)
set(CMAKE_CXX_EXTENSIONS    ON)
set(CMAKE_CXX_STANDARD      11)

include_directories(${PROJECT_SOURCE_DIR})

add_executable(${PROJECT_NAME} ${src} 
aec/delay_estimator/delay_estimator_wrapper.c 
aec/delay_estimator/delay_estimator.c 
aec/adp_filter_coeff_factor.c 
aec/aec_control.cc 
aec/aec_core.c 
subband/analy_synth/filterbank.c 
utility/fft/bitrv2.c 
utility/fft/bitrv2conj.c 
utility/fft/cdft_8g.c 
utility/fft/fft_wrapper.c 
utility/fft/ft_ops.c 
utility/fft/rdft_8g_init.c 
utility/fft/rdft_8g_make_table.c 
utility/fft/rdft_8g.c 
utility/ringbuffer/ring_buffer.c 
app/test1.cc)
target_link_libraries(${PROJECT_NAME} pthread)
