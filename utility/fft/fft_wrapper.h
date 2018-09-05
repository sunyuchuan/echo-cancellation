#ifndef _FFT_WRAPPER_H_
#define _FFT_WRAPPER_H_

#include "utility/fft/fft_config.h"

int rdft_create(FFT_Config** fft_conf, short fft_len);
int rdft_free(FFT_Config* fft_conf);

#endif
