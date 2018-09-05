#ifndef _FFT_CONFIGURE_H_
#define _FFT_CONFIGURE_H_

typedef struct _FFT_Config {
    int *ip;
    float *w;
    int nw;
    int nc;
    int fftLen;
} FFT_Config;

#endif