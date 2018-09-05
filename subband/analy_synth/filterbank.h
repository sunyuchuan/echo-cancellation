#ifndef _FILTER_BANK_H_
#define _FILTER_BANK_H_

#include "subband/analy_synth/filterbank_control.h"
#include "utility/fft/fft_config.h"

int DftFilterBankCreate(FilterBankControl **fb_inst);
int DftFilterBankInit(FilterBankControl *fb_inst);
int DftFilterBankAnalysis(FilterBankControl *fb_inst, FFT_Config *fft_config,
                          float *in_buf, short sample_num);
int DftFilterBankSynthesis(FilterBankControl *fb_inst, FFT_Config *fft_config,
                           float *out_buf, short *sample_num);
int DftFilterBankReset(FilterBankControl *fb_inst);
int DftFilterBankFree(FilterBankControl *fb_inst);

#endif