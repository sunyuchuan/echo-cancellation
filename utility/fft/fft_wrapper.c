#include "utility/fft/fft_wrapper.h"
#include <stdlib.h>
#include <string.h>
#include "aec/aec_defines.h"

int rdft_create(FFT_Config** fft_conf, short fft_len) {
    FFT_Config* fft_conf_t;

    fft_conf_t = (FFT_Config*)malloc(sizeof(FFT_Config));
    if (fft_conf_t == NULL) {
        return -1;
    }
    fft_conf_t->fftLen = fft_len;
    fft_conf_t->ip = (int*)malloc(sizeof(int) * (fft_len >> 2) * 5);
    if (fft_conf_t->ip == NULL) {
        return -1;
    }
    memset(fft_conf_t->ip, 0, sizeof(int) * (fft_len >> 2) * 5);
    fft_conf_t->w = (float*)malloc(sizeof(float) * (fft_len + 1));
    if (fft_conf_t->w == NULL) {
        return -1;
    }
    memset(fft_conf_t->w, 0, sizeof(float) * (fft_len + 1));
    fft_conf_t->nc = 0;
    fft_conf_t->nw = 0;

    *fft_conf = fft_conf_t;

    return 0;
}

int rdft_free(FFT_Config* fft_conf) {
    if (fft_conf == NULL) {
        return -1;
    }

    if (fft_conf->w != NULL) {
        free(fft_conf->w);
    }

    if (fft_conf->ip != NULL) {
        free(fft_conf->ip);
    }

    free(fft_conf);
    return 0;
}
