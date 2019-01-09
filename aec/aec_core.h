#ifndef _AEC_CORE_H_
#define _AEC_CORE_H_

#include "subband/analy_synth/filterbank_control.h"
#include "utility/fft/fft_config.h"

void CalcAbsValue(FilterBankControl *fb_inst, float *inv_abs_val,
                  float *abs_val);

void UpdateFarHistory(int *pos, float *abs_hist, float *abs_spectrum,
                      float *sepc_hist, float *spec);

float *AlignedFarend(int pos, float *far_hist, int delay);

void SetNonlinearGain(float curr_level, float *gain, float min_level);

int AecDeecho(float *far_signal, float *near_signal, float *far_frame,
              float *near_frame, float *filter, float *echo, float *error,
              float *abs_near, float *Rss, float *Rdd, float *Ree);

int AecResidualEchoCancellation(float *input_res, float *input_echo,
                                FFT_Config *fft_conf, float *pow_res,
                                float *pow_echo, float *crs_pow_rd,
                                float *adpt_pow_echo,
                                float *prev_enchanced_sqrd, float *enhanced,
                                short *first_frame, float nonlinear_gain);

int AecResidualEchoNN(float *input_res, float *input_echo, float *input_aligned_far,
		  float *concat_buf,float *nn_layer_buf, FFT_Config *fft_conf,
		  float *log_table,int *exp_table,int exp_precision, float *enhanced);

#endif
