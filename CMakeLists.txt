cmake_minimum_required(VERSION 3.4.1)

project(aec)

set (aec_sources
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
        )

add_library(${PROJECT_NAME} STATIC ${aec_sources})

set(CMAKE_C_EXTENSIONS      ON)
set(CMAKE_C_STANDARD        99)
set(CMAKE_CXX_EXTENSIONS    ON)
set(CMAKE_CXX_STANDARD      11)
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} ${CMAKE_CXX_FLAGS} -O3")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
