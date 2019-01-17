#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "subband/analy_synth/filterbank_control.h"
#include "aec/aec_defines.h"
#include "utility/fft/fft_config.h"
#include "utility/fft/rdft_8g.h"
#include "utility/math/fast_math.h"
#include "aec/post_process/sqrt_hanning_win.h"
#include "aec/post_process/aec_post_net.h"
#include "aec/post_process/hanning_win.h"
#include "aec/post_process/post_process_config.h"
#include "aec/aec_core.h"


#define ALPHA (-0.5f)
#define BETA (0.01f)
#define ALPHA_PLUS_1 (ALPHA + 1.0f)
#define ITERM1 ((1.0f - ALPHA) / (2.0f * 2.0f))  //(1-alpha)/(2*Adp_filter_length)

// post process
#define LAMDA (0.6f)
#define GAMMA (0.95f)

#define LOW_FREQUENCY_LEN 880

void CalcAbsValue(FilterBankControl *fb_inst, float *inv_abs_val,float *abs_val)
{
    short j, k;
    float *data1, tmp1, tmp2;

    data1 = &(fb_inst->analysis_fft_buf[0]);

    tmp1 = fabs(data1[0]);
    tmp2 = _reciprocal_sqrt(tmp1 + 1e-10f);
    abs_val[0] = tmp1;
    inv_abs_val[0] = tmp2 * tmp2;

    tmp1 = fabs(data1[1]);
    tmp2 = _reciprocal_sqrt(tmp1 + 1e-10f);
    abs_val[1] = tmp1;
    inv_abs_val[1] = tmp2 * tmp2;
    for (j = 2, k = 2; j < SUBBAND_NUM; j += 2, k++) {
        tmp1 = data1[j] * data1[j];
        tmp2 = tmp1 + data1[j + 1] * data1[j + 1];
        inv_abs_val[k] = _reciprocal_sqrt(tmp2 + 1e-10f);
        abs_val[k] = tmp2 * inv_abs_val[k];
    }
}

void UpdateFarHistory(int *pos, float *abs_hist, float *abs_spectrum,
                      float *sepc_hist, float *spec) {
    (*pos)++;
    if (*pos >= MAX_DELAY) {
        *pos = 0;
    }
    // Update far end spectrum buffer
    memcpy(&(abs_hist[*pos * PART_LEN1]), abs_spectrum,
           sizeof(float) * PART_LEN1);
    memcpy(&(sepc_hist[*pos * SUBBAND_NUM]), spec, sizeof(float) * SUBBAND_NUM);
}

float *AlignedFarend(int pos, float *spec_hist, int delay) {
    int buffer_position = pos - delay;

    // Check buffer position
    if (buffer_position < 0) {
        buffer_position += MAX_DELAY;
    }
    // Return far end spectrum
    return &(spec_hist[buffer_position * SUBBAND_NUM]);
}

void SetNonlinearGain(float curr_level, float *gain, float min_level) {
    if (curr_level <= min_level) {
        *gain = 1.0f;
    } else if (curr_level < (min_level + 0.1f)) {
        *gain = 0.8f;
    } else if (curr_level < (min_level + 0.2f)) {
        *gain = 0.7f;
    } else if (curr_level < (min_level + 0.3f)) {
        *gain = 0.6f;
    }else{
	*gain = 0.5f;
    }
}

int AecDeecho(float *far_signal, float *near_signal, float *far_frame,
              float *near_frame, float *filter, float *echo, float *error,
              float *abs_near, float *Rss, float *Rdd, float *Ree) {
    short i, j, k, subband_num;
    float tmpno1, tmpno2, tmpno3, tmpno4, tmpno5, tmpno6;
    float mu, fsum, k_l1, k_l2, k_l3, delta, kxr1, kxi1, kxr2, kxi2, kxr3, kxi3,
        cr, ci;

    error[0] = near_signal[0];
    error[1] = near_signal[1];
    echo[0] = 0.0f;
    echo[1] = 0.0f;
#if AEC_POST_PROCESSING_NN
	subband_num = 92;
#else
	subband_num = SUBBAND_NUM;
#endif

    for (i = 2, j = 6, k = 2; i < subband_num; i += 2, j += 6, k++) {
        // update far frame
        far_frame[j] = far_frame[j + 2];
        far_frame[j + 1] = far_frame[j + 3];
        far_frame[j + 2] = far_frame[j + 4];
        far_frame[j + 3] = far_frame[j + 5];
        far_frame[j + 4] = far_signal[i];
        far_frame[j + 5] = -far_signal[i + 1];
        // estimated echo
        tmpno1 = far_frame[j] * filter[j];
        tmpno2 = far_frame[j + 2] * filter[j + 2];
        tmpno3 = far_frame[j + 4] * filter[j + 4] + tmpno1;
        tmpno4 = far_frame[j + 1] * filter[j + 1] + tmpno2;
        tmpno5 = far_frame[j + 3] * filter[j + 3] + tmpno3;
        tmpno6 = far_frame[j + 5] * filter[j + 5] + tmpno4;
        echo[i] = tmpno5 + tmpno6;  // real

        tmpno1 = far_frame[j] * filter[j + 1];
        tmpno2 = far_frame[j + 1] * filter[j];  //-
        tmpno3 = far_frame[j + 2] * filter[j + 3] + tmpno1;
        tmpno4 = far_frame[j + 3] * filter[j + 2] + tmpno2;  //-
        tmpno5 = far_frame[j + 4] * filter[j + 5] + tmpno3;
        tmpno6 = far_frame[j + 5] * filter[j + 4] + tmpno4;  //-
        echo[i + 1] = tmpno5 - tmpno6;                       // imag
        // residual signal
        error[i] = near_signal[i] - echo[i];
        error[i + 1] = near_signal[i + 1] - echo[i + 1];
        // spectrum smooth
        tmpno2 = abs_near[k] * abs_near[k];
        tmpno3 = echo[i] * echo[i] + echo[i + 1] * echo[i + 1];
        tmpno5 = error[i] * error[i] + error[i + 1] * error[i + 1];
        tmpno1 = Rss[k] - tmpno2;
        tmpno4 = Rdd[k] - tmpno3;
        tmpno6 = Ree[k] - tmpno5;
        Rss[k] = 0.6f * tmpno1 + tmpno2;  // 0.9f->0.6f,20170329
        Rdd[k] = 0.6f * tmpno4 + tmpno3;
        Ree[k] = 0.6f * tmpno6 + tmpno5;
        // mu calculation
        tmpno3 = fabs(Rss[k] - Rdd[k]);
        tmpno2 = _reciprocal_sqrt(
            Ree[k] + 0.005f);  // 1/sqrt(Ree),0.001f->0.01f,20170110
                               // 0.01f->0.0005f,20170329
        tmpno4 = _reciprocal_sqrt(tmpno3 + 1e-10f);  // 1/sqrt(Rss - Rdd)
        tmpno5 = tmpno4 * tmpno3;                    // sqrt(Rss - Rdd)
        tmpno6 = fabs(1.0f - tmpno5 * tmpno2);       // sqrt((Rss - Rdd)/Ree)
        mu = fmin(tmpno6, 1.0f);
        // tendency matrix
        tmpno1 = filter[j] * filter[j];
        tmpno2 = filter[j + 1] * filter[j + 1] + tmpno1;
        tmpno3 = filter[j + 2] * filter[j + 2];
        tmpno4 = filter[j + 3] * filter[j + 3] + tmpno3;
        tmpno5 = filter[j + 4] * filter[j + 4];
        tmpno6 = filter[j + 5] * filter[j + 5] + tmpno5;
        // tmpno2 += 1e-30f;
        fsum = (tmpno2 + tmpno4 + tmpno6) * 2.0f;
        tmpno5 = _reciprocal(fsum + 1e-10f);
        tmpno3 = tmpno5 * ALPHA_PLUS_1;
        tmpno1 = 1.0f - mu;
        k_l1 = tmpno3 * tmpno2 + ITERM1;
        tmpno5 = tmpno1 * tmpno1;
        k_l2 = tmpno3 * tmpno4 + ITERM1;
        delta = tmpno5 * BETA;
        k_l3 = tmpno3 * tmpno6 + ITERM1;
        // increment of W numerator
        kxr1 = k_l1 * far_frame[j];
        kxi1 = k_l1 * far_frame[j + 1];
        kxr2 = k_l2 * far_frame[j + 2];
        kxi2 = k_l2 * far_frame[j + 3];
        kxr3 = k_l3 * far_frame[j + 4];
        kxi3 = k_l3 * far_frame[j + 5];
        // inverse regularization
        tmpno1 = kxr1 * far_frame[j] + delta;
        tmpno2 = kxi1 * far_frame[j + 1];
        tmpno3 = kxr2 * far_frame[j + 2] + tmpno1;
        tmpno4 = kxi2 * far_frame[j + 3] + tmpno2;
        tmpno5 = kxr3 * far_frame[j + 4] + tmpno3;
        tmpno6 = kxi3 * far_frame[j + 5] + tmpno4;
        tmpno1 = mu * error[i];
        tmpno3 = _reciprocal(tmpno5 + tmpno6 + 1e-10f);
        tmpno2 = mu * error[i + 1];
        cr = tmpno1 * tmpno3;
        ci = tmpno2 * tmpno3;

        tmpno1 = kxr1 * cr;
        tmpno2 = kxi1 * ci;
        tmpno4 = kxr1 * ci;
        tmpno3 = tmpno1 - tmpno2;
        tmpno5 = kxi1 * cr;
        filter[j] += tmpno3;
        tmpno1 = tmpno4 + tmpno5;
        tmpno2 = kxr2 * cr;
        tmpno3 = kxi2 * ci;
        filter[j + 1] += tmpno1;
        tmpno4 = kxr2 * ci;
        tmpno5 = tmpno2 - tmpno3;
        tmpno1 = kxi2 * cr;
        filter[j + 2] += tmpno5;
        tmpno2 = tmpno4 + tmpno1;
        filter[j + 3] += tmpno2;
        tmpno3 = kxr3 * cr;
        tmpno4 = kxi3 * ci;
        tmpno5 = kxr3 * ci;
        tmpno1 = tmpno3 - tmpno4;
        tmpno6 = kxi3 * cr;
        filter[j + 4] += tmpno1;
        tmpno2 = tmpno5 + tmpno6;
        filter[j + 5] += tmpno2;
    }
#if AEC_POST_PROCESSING_NN
	memcpy(&error[subband_num],&near_signal[subband_num],sizeof(float)*(SUBBAND_NUM-subband_num));
	memset(&echo[subband_num],0,sizeof(float)*(SUBBAND_NUM-subband_num));
#endif

    return 0;
}

int AecResidualEchoCancellation(float *input_res, float *input_echo,
                                FFT_Config *fft_conf, float *pow_res,
                                float *pow_echo, float *crs_pow_rd,
                                float *adpt_pow_echo,
                                float *prev_enchanced_sqrd, float *enhanced,
                                short *first_frame, float nonlinear_gain) {
    short i, j, k, p;
    float *in_res = input_res, *in_echo = input_echo;
    float tmpno1, tmpno2, tmpno3, tmpno4, tmpno5, tmpno6, tmpno7, tmpno8,
        invPbb, eta, ksi_dd, gain, res_fft_buf[POST_FFT_LEN],
        echo_fft_buf[POST_FFT_LEN];
    // update input data and window frame
    for (i = 0; i < POST_FFT_LEN; i += 2) {
        res_fft_buf[i] = sqrt_hanning_win[i] * in_res[i];
        res_fft_buf[i + 1] = sqrt_hanning_win[i + 1] * in_res[i + 1];
        echo_fft_buf[i] = sqrt_hanning_win[i] * in_echo[i];
        echo_fft_buf[i + 1] = sqrt_hanning_win[i + 1] * in_echo[i + 1];
    }
    // process data

    rdft(POST_FFT_LEN, &(fft_conf->nc), &(fft_conf->nw), 1, res_fft_buf,
         fft_conf->ip, fft_conf->w);
    rdft(POST_FFT_LEN, &(fft_conf->nc), &(fft_conf->nw), 1, echo_fft_buf,
         fft_conf->ip, fft_conf->w);

    memset(&res_fft_buf[LOW_FREQUENCY_LEN], 0,
           sizeof(float) * (POST_FFT_LEN - LOW_FREQUENCY_LEN));

    if (*first_frame == 1) {
        for (j = 2, k = 2; j < POST_FFT_LEN; j += 2, k++) {
            tmpno1 = res_fft_buf[j] * res_fft_buf[j];
            tmpno2 = echo_fft_buf[j] * echo_fft_buf[j];
            tmpno3 = res_fft_buf[j] * echo_fft_buf[j];
            tmpno4 = -res_fft_buf[j] * echo_fft_buf[j + 1];
            pow_res[k] = tmpno1 + res_fft_buf[j + 1] * res_fft_buf[j + 1];
            pow_echo[k] = tmpno2 + echo_fft_buf[j + 1] * echo_fft_buf[j + 1];
            adpt_pow_echo[k] = pow_echo[k];
            crs_pow_rd[j] = res_fft_buf[j + 1] * echo_fft_buf[j + 1] + tmpno3;
            crs_pow_rd[j + 1] = res_fft_buf[j + 1] * echo_fft_buf[j] + tmpno4;
            tmpno1 = crs_pow_rd[j] * crs_pow_rd[j];
            tmpno2 = crs_pow_rd[j + 1] * crs_pow_rd[j + 1] + tmpno1;  //|Prd|^2
            // tmpno3 = _reciprocal(tmpno2*adpt_pow_echo[k]+1e-12f);
            // tmpno4 = pow_res[k]*pow_echo[k];
            tmpno3 = _reciprocal(pow_res[k] * pow_echo[k] + 1e-12f);
            tmpno4 = tmpno2 * adpt_pow_echo[k];
            invPbb = _reciprocal(tmpno3 * tmpno4 + 1e-12f);
            eta = fmax(pow_res[k] * invPbb, 1.0f);
            ksi_dd = fmax(eta - 1.0f, 0.0001f);
            tmpno1 = _reciprocal_sqrt(ksi_dd + 1.0f);
            tmpno2 = _reciprocal_sqrt(ksi_dd + 1e-12f);
            tmpno3 = ksi_dd * tmpno2;
            gain = tmpno3 * tmpno1;
            // gain = 1.0f;
            enhanced[j] = res_fft_buf[j] * gain * 0.75f;
            enhanced[j + 1] = res_fft_buf[j + 1] * gain * 0.75f;
            prev_enchanced_sqrd[k] =
                enhanced[j] * enhanced[j] + enhanced[j + 1] * enhanced[j + 1];
        }
        (*first_frame)++;
    } else if (*first_frame < 88)  // one second for delay estimation
                                   // convergence, make the gain more aggressive
    {
        for (j = 2, k = 2; j < POST_FFT_LEN; j += 2, k++) {
            tmpno1 = res_fft_buf[j] * res_fft_buf[j];
            tmpno2 = echo_fft_buf[j] * echo_fft_buf[j];
            tmpno3 = res_fft_buf[j] * echo_fft_buf[j];
            tmpno4 = -res_fft_buf[j] * echo_fft_buf[j + 1];
            tmpno5 = res_fft_buf[j + 1] * res_fft_buf[j + 1] + tmpno1;
            tmpno6 = echo_fft_buf[j + 1] * echo_fft_buf[j + 1] + tmpno2;
            tmpno1 = res_fft_buf[j + 1] * echo_fft_buf[j + 1] + tmpno3;
            tmpno2 = res_fft_buf[j + 1] * echo_fft_buf[j] + tmpno4;

            tmpno3 = pow_res[k] - tmpno5;
            tmpno4 = pow_echo[k] - tmpno6;
            tmpno7 = crs_pow_rd[j] - tmpno1;
            tmpno8 = crs_pow_rd[j + 1] - tmpno2;
            pow_res[k] = tmpno3 * LAMDA + tmpno5;
            pow_echo[k] = tmpno4 * LAMDA + tmpno6;
            crs_pow_rd[j] = tmpno7 * LAMDA + tmpno1;
            crs_pow_rd[j + 1] = tmpno8 * LAMDA + tmpno2;
            tmpno3 = _reciprocal(pow_echo[k] + 1e-12f);
            tmpno4 = fmin(tmpno3 * tmpno6, 1.0f);
            tmpno5 = tmpno6 - adpt_pow_echo[k];
            adpt_pow_echo[k] += tmpno4 * tmpno5;
            tmpno1 = crs_pow_rd[j] * crs_pow_rd[j];
            tmpno2 = crs_pow_rd[j + 1] * crs_pow_rd[j + 1] + tmpno1;  //|Prd|^2
            /*tmpno3 = _reciprocal(tmpno2*adpt_pow_echo[k]+1e-12f);
            tmpno4 = pow_res[k]*pow_echo[k];*/

            // Crd^2*Pee -> Crd*Pee,20170110
            tmpno3 = _reciprocal(pow_res[k] * pow_echo[k] + 1e-12f);
            tmpno4 = tmpno2 * adpt_pow_echo[k];
            invPbb = _reciprocal(tmpno3 * tmpno4 + 1e-12f);

            //			tmpno3 = adpt_pow_echo[k]*adpt_pow_echo[k];
            //			tmpno5 = pow_res[k]*pow_echo[k];
            //			tmpno4 = _reciprocal_sqrt(tmpno2*tmpno3+1e-12f);
            //			tmpno6 = _reciprocal_sqrt(tmpno5+1e-12f);
            //			tmpno3 = tmpno5*tmpno6;
            //			invPbb = tmpno3*tmpno4;
            // 20170110, modification end

            tmpno7 = fmax(pow_res[k] * invPbb, 1.0001f);
            eta = fmin(tmpno7, 10000.0f);  // 10lg(10000)=40
            tmpno5 = eta - 1.0f;
            tmpno6 = prev_enchanced_sqrd[k] * invPbb - tmpno5;
            ksi_dd = GAMMA * tmpno6 + tmpno5;
            tmpno1 = 1.0f + ksi_dd;
            tmpno2 = ksi_dd * ksi_dd;
            ;
            tmpno3 = tmpno2 * eta;
            tmpno4 = tmpno1 * tmpno1;
            tmpno5 = _reciprocal_sqrt(tmpno3 + 1e-12f);
            tmpno6 = _reciprocal_sqrt(tmpno3 + tmpno4 + 1e-12f);
            tmpno1 = tmpno5 * tmpno3;
            gain = tmpno1 * tmpno6;
            enhanced[j] = res_fft_buf[j] * gain * 0.75f;
            enhanced[j + 1] = res_fft_buf[j + 1] * gain * 0.75f;
            prev_enchanced_sqrd[k] =
                enhanced[j] * enhanced[j] + enhanced[j + 1] * enhanced[j + 1];
        }
        (*first_frame)++;
    } else {
        for (j = 2, k = 2; j < POST_FFT_LEN; j += 2, k++) {
            tmpno1 = res_fft_buf[j] * res_fft_buf[j];
            tmpno2 = echo_fft_buf[j] * echo_fft_buf[j];
            tmpno3 = res_fft_buf[j] * echo_fft_buf[j];
            tmpno4 = -res_fft_buf[j] * echo_fft_buf[j + 1];
            tmpno5 = res_fft_buf[j + 1] * res_fft_buf[j + 1] + tmpno1;
            tmpno6 = echo_fft_buf[j + 1] * echo_fft_buf[j + 1] + tmpno2;
            tmpno1 = res_fft_buf[j + 1] * echo_fft_buf[j + 1] + tmpno3;
            tmpno2 = res_fft_buf[j + 1] * echo_fft_buf[j] + tmpno4;

            tmpno3 = pow_res[k] - tmpno5;
            tmpno4 = pow_echo[k] - tmpno6;
            tmpno7 = crs_pow_rd[j] - tmpno1;
            tmpno8 = crs_pow_rd[j + 1] - tmpno2;
            pow_res[k] = tmpno3 * LAMDA + tmpno5;
            pow_echo[k] = tmpno4 * LAMDA + tmpno6;
            crs_pow_rd[j] = tmpno7 * LAMDA + tmpno1;
            crs_pow_rd[j + 1] = tmpno8 * LAMDA + tmpno2;
            tmpno3 = _reciprocal(pow_echo[k] + 1e-12f);
            tmpno4 = fmin(tmpno3 * tmpno6, 1.0f);
            tmpno5 = tmpno6 - adpt_pow_echo[k];
            adpt_pow_echo[k] += tmpno4 * tmpno5;
            tmpno1 = crs_pow_rd[j] * crs_pow_rd[j];
            tmpno2 = crs_pow_rd[j + 1] * crs_pow_rd[j + 1] + tmpno1;  //|Prd|^2
            /*tmpno3 = _reciprocal(tmpno2*adpt_pow_echo[k]+1e-12f);
            tmpno4 = pow_res[k]*pow_echo[k];*/

            // Crd^2*Pee -> Crd*Pee,20170110
            tmpno3 = _reciprocal(pow_res[k] * pow_echo[k] + 1e-12f);
            tmpno4 = tmpno2 * adpt_pow_echo[k];
            invPbb = _reciprocal(tmpno3 * tmpno4 + 1e-12f);

            //			tmpno3 = adpt_pow_echo[k]*adpt_pow_echo[k];
            //			tmpno5 = pow_res[k]*pow_echo[k];
            //			tmpno4 = _reciprocal_sqrt(tmpno2*tmpno3+1e-12f);
            //			tmpno6 = _reciprocal_sqrt(tmpno5+1e-12f);
            //			tmpno3 = tmpno5*tmpno6;
            //			invPbb = tmpno3*tmpno4;
            // 20170110, modification end

            tmpno7 = fmax(pow_res[k] * invPbb, 1.0001f);
            eta = fmin(tmpno7, 10000.0f);  // 10lg(10000)=40
            tmpno5 = eta - 1.0f;
            tmpno6 = prev_enchanced_sqrd[k] * invPbb - tmpno5;
            ksi_dd = GAMMA * tmpno6 + tmpno5;
            tmpno1 = 1.0f + ksi_dd;
            tmpno2 = ksi_dd * ksi_dd;
            ;
            tmpno3 = tmpno2 * eta;
            tmpno4 = tmpno1 * tmpno1;
            tmpno5 = _reciprocal_sqrt(tmpno3 + 1e-12f);
            tmpno6 = _reciprocal_sqrt(tmpno3 + tmpno4 + 1e-12f);
            tmpno1 = tmpno5 * tmpno3;
            gain = tmpno1 * tmpno6;
            enhanced[j] = res_fft_buf[j] * gain * nonlinear_gain;
            enhanced[j + 1] = res_fft_buf[j + 1] * gain * nonlinear_gain;
            prev_enchanced_sqrd[k] =
                enhanced[j] * enhanced[j] + enhanced[j + 1] * enhanced[j + 1];
        }
    }
    enhanced[0] = 0.0f;
    enhanced[1] = 0.0f;

    rdft(POST_FFT_LEN, &(fft_conf->nc), &(fft_conf->nw), -1, enhanced,
         fft_conf->ip, fft_conf->w);
    for (p = 0; p < POST_FFT_LEN; p += 4) {
        enhanced[p] *= sqrt_hanning_win[p];
        enhanced[p + 1] *= sqrt_hanning_win[p + 1];
        enhanced[p + 2] *= sqrt_hanning_win[p + 2];
        enhanced[p + 3] *= sqrt_hanning_win[p + 3];
    }

    return 0;
}

static void AecPostProcess_DenseLayer_555X256_ActivationTanh(const float input[555],const float output[256], unsigned int *exp_table, int exp_precision)
{
	int i,j;
	float sum1 = 0.0f,sum2 = 0.0f,sum3 = 0.0f,sum4 = 0.0f,sum5 = 0.0f,sum6 = 0.0f,sum7 = 0.0f,sum8 = 0.0f;
	const float *w1,*w2,*w3,*w4,*w5,*w6,*w7,*w8,*in1,*in2,*in3,*in4,*in5,*in6,*in7,*in8;

	in1 = &input[0];
	in2 = &input[69];
	in3 = &input[138];
	in4 = &input[207];
	in5 = &input[276];
	in6 = &input[345];
	in7 = &input[414];
	in8 = &input[483];
	for(i=0;i<256;i+=1)
	{
		w1 = &W_DenseLayer1[i][0];
		w2 = &W_DenseLayer1[i][69];
		w3 = &W_DenseLayer1[i][138];
		w4 = &W_DenseLayer1[i][207];
		w5 = &W_DenseLayer1[i][276];
		w6 = &W_DenseLayer1[i][345];
		w7 = &W_DenseLayer1[i][414];
		w8 = &W_DenseLayer1[i][483];
		for(j=0;j<69;j+=1)
		{
			sum1+=in1[j] * w1[j];
			sum2+=in2[j] * w2[j];
			sum3+=in3[j] * w3[j];
			sum4+=in4[j] * w4[j];
			sum5+=in5[j] * w5[j];
			sum6+=in6[j] * w6[j];
			sum7+=in7[j] * w7[j];
			sum8+=in8[j] * w8[j];
		}
		sum1+=input[552]*w1[552];
		sum2+=input[553]*w1[553];
		sum3+=input[554]*w1[554];
		sum4+=W_DenseLayer1_Bias[i];
		sum5+=(sum1 + sum2 + sum3 + sum4 + sum6 + sum7 + sum8);
		output[i] = _tanh(sum5,exp_table,AEC_EXP_PRECISION);
		sum1 = 0.f;
		sum2 = 0.f;
		sum3 = 0.f;
		sum4 = 0.f;
		sum5 = 0.f;
		sum6 = 0.f;
		sum7 = 0.f;
		sum8 = 0.f;
	}
}

static void AecPostProcess_DenseLayer_256X192_ActivationTanh(const float input[256],float output[192],unsigned int *exp_table, int exp_precision)
{
	int i,j;
	float sum1 = 0.0f,sum2 = 0.0f,sum3 = 0.0f,sum4 = 0.0f,sum5 = 0.0f,sum6 = 0.0f,sum7 = 0.0f,sum8 = 0.0f;
	const float *w1,*w2,*w3,*w4,*w5,*w6,*w7,*w8,*in1,*in2,*in3,*in4,*in5,*in6,*in7,*in8;

	in1 = &input[0];
	in2 = &input[32];
	in3 = &input[64];
	in4 = &input[96];
	in5 = &input[128];
	in6 = &input[160];
	in7 = &input[192];
	in8 = &input[224];
	for(i=0;i<192;i+=1)
	{
		w1 = &W_DenseLayer2[i][0];
		w2 = &W_DenseLayer2[i][32];
		w3 = &W_DenseLayer2[i][64];
		w4 = &W_DenseLayer2[i][96];
		w5 = &W_DenseLayer2[i][128];
		w6 = &W_DenseLayer2[i][160];
		w7 = &W_DenseLayer2[i][192];
		w8 = &W_DenseLayer2[i][224];
		for(j=0;j<32;j+=1)
		{
			sum1+=in1[j] * w1[j];
			sum2+=in2[j] * w2[j];
			sum3+=in3[j] * w3[j];
			sum4+=in4[j] * w4[j];
			sum5+=in5[j] * w5[j];
			sum6+=in6[j] * w6[j];
			sum7+=in7[j] * w7[j];
			sum8+=in8[j] * w8[j];
		}
		sum5+=(sum1 + sum2 + sum3 + sum4 + sum6 + sum7 + sum8 + W_DenseLayer2_Bias[i]);
		output[i] = _tanh(sum5,exp_table,AEC_EXP_PRECISION);
		sum1 = 0.f;
		sum2 = 0.f;
		sum3 = 0.f;
		sum4 = 0.f;
		sum5 = 0.f;
		sum6 = 0.f;
		sum7 = 0.f;
		sum8 = 0.f;
	}
}

static void AecPostProcess_DenseLayer_192X192_ActivationTanh(const float input[192],const float output[192], unsigned int *exp_table, int exp_precision)
{
	int i,j;
	float sum1 = 0.0f,sum2 = 0.0f,sum3 = 0.0f,sum4 = 0.0f,sum5 = 0.0f,sum6 = 0.0f,sum7 = 0.0f,sum8 = 0.0f;
	const float *w1,*w2,*w3,*w4,*w5,*w6,*w7,*w8,*in1,*in2,*in3,*in4,*in5,*in6,*in7,*in8;

	in1 = &input[0];
	in2 = &input[24];
	in3 = &input[48];
	in4 = &input[72];
	in5 = &input[96];
	in6 = &input[120];
	in7 = &input[144];
	in8 = &input[168];
	for(i=0;i<192;i+=1)
	{
		w1 = &W_DenseLayer3[i][0];
		w2 = &W_DenseLayer3[i][24];
		w3 = &W_DenseLayer3[i][48];
		w4 = &W_DenseLayer3[i][72];
		w5 = &W_DenseLayer3[i][96];
		w6 = &W_DenseLayer3[i][120];
		w7 = &W_DenseLayer3[i][144];
		w8 = &W_DenseLayer3[i][168];
		for(j=0;j<24;j+=1)
		{
			sum1+=in1[j] * w1[j];
			sum2+=in2[j] * w2[j];
			sum3+=in3[j] * w3[j];
			sum4+=in4[j] * w4[j];
			sum5+=in5[j] * w5[j];
			sum6+=in6[j] * w6[j];
			sum7+=in7[j] * w7[j];
			sum8+=in8[j] * w8[j];
		}
		sum5+=(sum1 + sum2 + sum3 + sum4 + sum6 + sum7 + sum8 + W_DenseLayer3_Bias[i]);
		output[i] = _tanh(sum5,exp_table,AEC_EXP_PRECISION);
		sum1 = 0.f;
		sum2 = 0.f;
		sum3 = 0.f;
		sum4 = 0.f;
		sum5 = 0.f;
		sum6 = 0.f;
		sum7 = 0.f;
		sum8 = 0.f;
	}
}

static void AecPostProcess_DenseLayer_192X185_ActivationTanh(const float input[192],const float output[185], unsigned int *exp_table, int exp_precision)
{
	int i,j;
	float sum1 = 0.0f,sum2 = 0.0f,sum3 = 0.0f,sum4 = 0.0f,sum5 = 0.0f,sum6 = 0.0f,sum7 = 0.0f,sum8 = 0.0f;
	const float *w1,*w2,*w3,*w4,*w5,*w6,*w7,*w8,*in1,*in2,*in3,*in4,*in5,*in6,*in7,*in8;

	in1 = &input[0];
	in2 = &input[24];
	in3 = &input[48];
	in4 = &input[72];
	in5 = &input[96];
	in6 = &input[120];
	in7 = &input[144];
	in8 = &input[168];
	for(i=0;i<185;i+=1)
	{
		w1 = &W_DenseLayer4[i][0];
		w2 = &W_DenseLayer4[i][24];
		w3 = &W_DenseLayer4[i][48];
		w4 = &W_DenseLayer4[i][72];
		w5 = &W_DenseLayer4[i][96];
		w6 = &W_DenseLayer4[i][120];
		w7 = &W_DenseLayer4[i][144];
		w8 = &W_DenseLayer4[i][168];
		for(j=0;j<24;j+=1)
		{
			sum1+=in1[j] * w1[j];
			sum2+=in2[j] * w2[j];
			sum3+=in3[j] * w3[j];
			sum4+=in4[j] * w4[j];
			sum5+=in5[j] * w5[j];
			sum6+=in6[j] * w6[j];
			sum7+=in7[j] * w7[j];
			sum8+=in8[j] * w8[j];
		}
		sum5+=(sum1 + sum2 + sum3 + sum4 + sum6 + sum7 + sum8 + W_DenseLayer4_Bias[i]);
		output[i] = _tanh(sum5,exp_table,AEC_EXP_PRECISION);
		sum1 = 0.f;
		sum2 = 0.f;
		sum3 = 0.f;
		sum4 = 0.f;
		sum5 = 0.f;
		sum6 = 0.f;
		sum7 = 0.f;
		sum8 = 0.f;
	}
}

int AecResidualEchoNN(float *input_res, float *input_echo, float *input_aligned_far,
						float *concat_buf,float *nn_layer_buf, FFT_Config *fft_conf,
						float *log_table,int *exp_table,int exp_precision, float *enhanced)
{
	short i,j,p,q,m,n,f1,f2,f3;
	float *in_res = input_res, *in_echo = input_echo, *in_far = input_aligned_far, *_3dim_buf,*max_nn_dim_buf, *output_buf = enhanced;
	float tmpno1,tmpno2,tmpno3,tmpno4,a1,a2,res_fft_buf[POST_FFT_LEN],echo_fft_buf[POST_FFT_LEN],far_fft_buf[POST_FFT_LEN];
	
	/*for(i=0;i<POST_FFT_LEN;i++)
	{
		f1 = (short)(input_res[i]*32767.0f);
		input_res[i] = ((float)f1)*(1.0f/32768.0f);
		f2 = (short)(input_echo[i]*32767.0f);
		input_echo[i] = ((float)f2)*(1.0f/32768.0f);
		f3 = (short)(input_aligned_far[i]*32767.0f);
		input_aligned_far[i] = ((float)f3)*(1.0f/32768.0f);
	}*/
	
	//update input data and window frame
	for(i=0;i<POST_FFT_LEN;i+=2)
	{
		//sum1+=(in_res[i]+in_res[i+1]);
		res_fft_buf[i] = hanning_win[i]*in_res[i];
		res_fft_buf[i+1] = hanning_win[i+1]*in_res[i+1];
		echo_fft_buf[i] = hanning_win[i]*in_echo[i];
		echo_fft_buf[i+1] = hanning_win[i+1]*in_echo[i+1];
		far_fft_buf[i] = hanning_win[i]*in_far[i];
		far_fft_buf[i+1] = hanning_win[i+1]*in_far[i+1];
		/*if(res_fft_buf[i]!=0.0f || res_fft_buf[i+1]!=0.0f)
		{
			sum2 = 0.0f;
		}*/
	}
	//process data
	rdft(POST_FFT_LEN,&(fft_conf->nc),&(fft_conf->nw),1,res_fft_buf,fft_conf->ip,fft_conf->w);
	rdft(POST_FFT_LEN,&(fft_conf->nc),&(fft_conf->nw),1,echo_fft_buf,fft_conf->ip,fft_conf->w);
	rdft(POST_FFT_LEN,&(fft_conf->nc),&(fft_conf->nw),1,far_fft_buf,fft_conf->ip,fft_conf->w);
	//concatenate data
	//input  = (0.5*ln(a*conj(a)) - mean)/global_max_val,mean = 4.939886,global_max_val = 1/0.08945976
	//		 = (0.5/global_max_val)*ln(a*conj(a)) - mean/global_max_val
	//		 = a1**ln(a*conj(a)) + a2
	a1 = 0.08945976f;
	a2 = -(-4.939886)*0.08945976f;
	_3dim_buf = concat_buf;
	max_nn_dim_buf = nn_layer_buf;
	//tmpno1 = logf(fabs(res_fft_buf[0])+1e-7f);
	tmpno1 = _ln(fabs(res_fft_buf[0])+1e-7f,log_table,AEC_LN_PRECISION);
	_3dim_buf[0] = tmpno1*a1 + a2;
	for(p=2,q=1;p<TARGET_DIM*2;p+=2,q++)
	{
		tmpno1 = res_fft_buf[p]*res_fft_buf[p];
		tmpno2 = res_fft_buf[p+1]*res_fft_buf[p+1] + tmpno1;
		tmpno3 = _reciprocal_sqrt_hp(tmpno2+1e-8f);
		//tmpno1 = logf(sqrt(tmpno2)+1e-7f);
		tmpno1 = _ln(tmpno3*tmpno2+1e-7f,log_table,AEC_LN_PRECISION);
		_3dim_buf[q] = tmpno1*a1 + a2;
	}
	//tmpno1 = logf(fabs(echo_fft_buf[0])+1e-7f);
	tmpno1 = _ln(fabs(echo_fft_buf[0])+1e-7f,log_table,AEC_LN_PRECISION);
	_3dim_buf[TARGET_DIM] = tmpno1*a1 + a2;
	for(i=2,j=1;i<TARGET_DIM*2;i+=2,j++)
	{
		tmpno1 = echo_fft_buf[i]*echo_fft_buf[i];
		tmpno2 = echo_fft_buf[i+1]*echo_fft_buf[i+1] + tmpno1;
		tmpno3 = _reciprocal_sqrt_hp(tmpno2+1e-8f);
		//tmpno1 = logf(sqrt(tmpno2)+1e-7f);
		tmpno1 = _ln(tmpno3*tmpno2+1e-7f,log_table,AEC_LN_PRECISION);
		_3dim_buf[TARGET_DIM+j] = tmpno1*a1 + a2;
	}
	//tmpno1 = logf(fabs(far_fft_buf[0])+1e-7f);
	tmpno1 = _ln(fabs(far_fft_buf[0])+1e-7f,log_table,AEC_LN_PRECISION);
	_3dim_buf[TARGET_DIM*2] = tmpno1*a1 + a2;
	for(m=2,n=1;m<TARGET_DIM*2;m+=2,n++)
	{
		tmpno1 = far_fft_buf[m]*far_fft_buf[m];
		tmpno2 = far_fft_buf[m+1]*far_fft_buf[m+1] + tmpno1;
		tmpno3 = _reciprocal_sqrt_hp(tmpno2+1e-8f);
		//tmpno1 = logf(sqrt(tmpno2)+1e-7f);
		tmpno1 = _ln(tmpno3*tmpno2+1e-7f,log_table,AEC_LN_PRECISION);
		_3dim_buf[2*TARGET_DIM+n] = tmpno1*a1 + a2;
	}

	AecPostProcess_DenseLayer_555X256_ActivationTanh(_3dim_buf,max_nn_dim_buf,exp_table,exp_precision);
	AecPostProcess_DenseLayer_256X192_ActivationTanh(max_nn_dim_buf,_3dim_buf,exp_table,exp_precision);
	AecPostProcess_DenseLayer_192X192_ActivationTanh(_3dim_buf,max_nn_dim_buf,exp_table,exp_precision);
	AecPostProcess_DenseLayer_192X185_ActivationTanh(max_nn_dim_buf,_3dim_buf,exp_table,exp_precision);

	output_buf[0] = res_fft_buf[0]*_3dim_buf[0];
	output_buf[1] = 0.0f;
	for(i=2,j=1;i<TARGET_DIM*2;i+=2,j++)
	{
		output_buf[i] = res_fft_buf[i]*_3dim_buf[j];
		output_buf[i+1] = res_fft_buf[i+1]*_3dim_buf[j];
	}

	memset(&output_buf[TARGET_DIM*2],0,(POST_FFT_LEN-TARGET_DIM*2)*sizeof(float));
	rdft(POST_FFT_LEN,&(fft_conf->nc),&(fft_conf->nw),-1,output_buf,fft_conf->ip,fft_conf->w);

	return 0;
}
