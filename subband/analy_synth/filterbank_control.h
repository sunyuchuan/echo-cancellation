#ifndef _FILTER_BANK_CONTROL_H_
#define _FILTER_BANK_CONTROL_H_

typedef struct _FilterBankControl {
    short subband_num;
    short oversample_rate;
    short prototype_filter_len;

    short analysis_buf_head_pos;
    float *analysis_buf;
    float *analysis_fft_buf;

    float *synthesis_fft_buf;
    float *synthesis_state_buf;
} FilterBankControl;

#endif