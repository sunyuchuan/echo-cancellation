#ifndef _AEC_DEFINES_H_
#define _AEC_DEFINES_H_

#include "aec/aec_config.h"

#ifdef SUBBAND_AEC_PROTOTYPE_FILTER_2048

#define SUBBAND_NUM 512  // 44k
#define OVERSAMPLE_RATE 2
#define PART_LEN (SUBBAND_NUM >> 1)
#define PART_LEN1 (PART_LEN + 1)  // even-stacking mode. Nyquist band included.
#define FAR_BUF_LEN (PART_LEN << 2)

#define FRAME_LEN SUBBAND_NUM
#define MAX_DELAY 50
#define ADPF_LEN (512)  // 44k
#define SUBBAND_FRAME_SHIFT PART_LEN
#define FAR_END_BUF_FRAME_NUM 50

#else

#define SUBBAND_NUM 256  // 44k
#define OVERSAMPLE_RATE 2
#define PART_LEN (SUBBAND_NUM >> 1)
#define PART_LEN1 (PART_LEN + 1)  // even-stacking mode. Nyquist band included.
#define FAR_BUF_LEN (PART_LEN << 2)

#define FRAME_LEN SUBBAND_NUM
#define MAX_DELAY 120
#define ADPF_LEN (384)  // 44k
#define SUBBAND_FRAME_SHIFT PART_LEN
#define FAR_END_BUF_FRAME_NUM 128

#endif

#define POST_FFT_LEN 1024

#endif
