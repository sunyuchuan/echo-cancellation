cmake_minimum_required(VERSION 2.8)

SET(CMAKE_BUILD_TYPE "Debug")  
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")  
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

project(echo_cancellation)

include_directories(${PROJECT_SOURCE_DIR})

add_executable(${PROJECT_NAME} ${src} 
aec/delay_estimator/delay_estimator_wrapper.c 
aec/delay_estimator/delay_estimator.c 
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
target_link_libraries(${PROJECT_NAME})
