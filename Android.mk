LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE    := EchoCancellation
LOCAL_SRC_FILES := aec/delay_estimator/delay_estimator_wrapper.c \
					aec/delay_estimator/delay_estimator.c \
					aec/adp_filter_coeff_factor.c \
					aec/aec_control.cc \
					aec/aec_core.c \
					subband/analy_synth/filterbank.c \
					utility/fft/bitrv2.c \
					utility/fft/bitrv2conj.c \
					utility/fft/cdft_8g.c \
					utility/fft/fft_wrapper.c \
					utility/fft/ft_ops.c \
					utility/fft/rdft_8g_init.c \
					utility/fft/rdft_8g_make_table.c \
					utility/fft/rdft_8g.c \
					utility/ringbuffer/ring_buffer.c
					
LOCAL_CFLAGS := -mfloat-abi=softfp -mfpu=neon -marm -ftree-vectorize -ffast-math -O3
LOCAL_CXXFLAGS := -std=c++11
LOCAL_DEFAULT_CPP_EXTENSION := '.c' '.cc'
LOCAL_LDLIBS := -llog -lm

include $(BUILD_STATIC_LIBRARY)
