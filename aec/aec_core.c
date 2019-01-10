#include "aec/aec_core.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aec/aec_defines.h"
#include "aec/post_process/sqrt_hanning_win.h"
#include "utility/fft/rdft_8g.h"
#include "utility/math/fast_math.h"
#include "aec/post_process/aec_post_net.h"
#include "aec/post_process/hanning_win.h"
#include "aec/post_process/post_process_config.h"

#define ALPHA (-0.5f)
#define BETA (0.01f)
#define ALPHA_PLUS_1 (ALPHA + 1.0f)
#define ITERM1 \
    ((1.0f - ALPHA) / (2.0f * 2.0f))  //(1-alpha)/(2*Adp_filter_length)

// post process
#define LAMDA (0.6f)
#define GAMMA (0.95f)

#define LOW_FREQUENCY_LEN 880

void CalcAbsValue(FilterBankControl *fb_inst, float *inv_abs_val,
                  float *abs_val) {
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

static void AecPostProcess_DenseLayer_555X512_ActivationTanh(float *input,float *output, int *exp_table, int exp_precision)
{
	int i,j;
	float sum1 = 0.0f,sum2 = 0.0f,sum3 = 0.0f,sum4 = 0.0f,sum5 = 0.0f,*w;

	for(i=0;i<512;i++)
	{
		w = &W_DenseLayer1[i][0];
		for(j=0;j<552;j+=4)
		{
			sum1+=input[j]*w[j];
			sum2+=input[j+1]*w[j+1];
			sum3+=input[j+2]*w[j+2];
			sum4+=input[j+3]*w[j+3];
			//sum5+=input[j+4]*w[j+4];
		}
		sum2+=input[552]*w[552];
		sum3+=input[553]*w[553];
		sum4+=input[554]*w[554];
		sum1+=(sum2 + sum3 + sum4 + sum5 + W_DenseLayer1_Bias[i]);
		output[i] = _tanh(sum1,exp_table,exp_precision);
		//output[i] = tanh(sum1);
		sum1 = 0.f;
		sum2 = 0.f;
		sum3 = 0.f;
		sum4 = 0.f;
		sum5 = 0.f;
	}
}

static void AecPostProcess_DenseLayer_512X256_ActivationTanh(float *input,float *output, int *exp_table, int exp_precision)
{
	int i,j;
	float sum1 = 0.0f,sum2 = 0.0f,sum3 = 0.0f,sum4 = 0.0f,*w;

	for(i=0;i<256;i++)
	{
		w = &W_DenseLayer2[i][0];
		for(j=0;j<512;j+=4)
		{
			sum1+=input[j]*w[j];
			sum2+=input[j+1]*w[j+1];
			sum3+=input[j+2]*w[j+2];
			sum4+=input[j+3]*w[j+3];
		}

		sum1+=(sum2 + sum3 + sum4 + W_DenseLayer2_Bias[i]);
		output[i] = _tanh(sum1,exp_table,exp_precision);
		//output[i] = tanh(sum1);
		sum1 = 0.f;
		sum2 = 0.f;
		sum3 = 0.f;
		sum4 = 0.f;
	}
}

static void AecPostProcess_DenseLayer_256X256_ActivationTanh(float *input,float *output, int *exp_table, int exp_precision)
{
	int i,j;
	float sum1 = 0.0f,sum2 = 0.0f,sum3 = 0.0f,sum4 = 0.0f,*w;

	for(i=0;i<256;i++)
	{
		w = &W_DenseLayer3[i][0];
		for(j=0;j<256;j+=4)
		{
			sum1+=input[j]*w[j];
			sum2+=input[j+1]*w[j+1];
			sum3+=input[j+2]*w[j+2];
			sum4+=input[j+3]*w[j+3];
		}
		sum1+=(sum2 + sum3 + sum4 + W_DenseLayer3_Bias[i]);
		output[i] = _tanh(sum1,exp_table,exp_precision);
		//output[i] = tanh(sum1);
		sum1 = 0.f;
		sum2 = 0.f;
		sum3 = 0.f;
		sum4 = 0.f;
	}

}

static void AecPostProcess_DenseLayer_256X185_ActivationTanh(float *input,float *output, int *exp_table, int exp_precision)
{
	int i,j;
	float sum1 = 0.0f,sum2 = 0.0f,sum3 = 0.0f,sum4 = 0.0f,*w;

	for(i=0;i<185;i++)
	{
		w = &W_DenseLayer4[i][0];
		for(j=0;j<256;j+=4)
		{
			sum1+=input[j]*w[j];
			sum2+=input[j+1]*w[j+1];
			sum3+=input[j+2]*w[j+2];
			sum4+=input[j+3]*w[j+3];
		}
		sum1+=(sum2 + sum3 + sum4 + W_DenseLayer4_Bias[i]);
		output[i] = _tanh(sum1,exp_table,exp_precision);
		//output[i] = tanh(sum1);
		sum1 = 0.f;
		sum2 = 0.f;
		sum3 = 0.f;
		sum4 = 0.f;
	}
}

int AecResidualEchoNN(float *input_res, float *input_echo, float *input_aligned_far,
						float *concat_buf,float *nn_layer_buf, FFT_Config *fft_conf,
						float *log_table,int *exp_table,int exp_precision, float *enhanced)
{
	short i,j,p,q,m,n,f1,f2,f3;
	float *in_res = input_res, *in_echo = input_echo, *in_far = input_aligned_far;
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
	//tmpno1 = logf(fabs(res_fft_buf[0])+1e-7f);
	tmpno1 = _ln(fabs(res_fft_buf[0])+1e-7f,log_table,AEC_LN_PRECISION);
	concat_buf[0] = tmpno1*a1 + a2;
	for(p=2,q=1;p<TARGET_DIM*2;p+=2,q++)
	{
		tmpno1 = res_fft_buf[p]*res_fft_buf[p];
		tmpno2 = res_fft_buf[p+1]*res_fft_buf[p+1] + tmpno1;
		tmpno3 = _reciprocal_sqrt_hp(tmpno2+1e-8f);
		//tmpno1 = logf(sqrt(tmpno2)+1e-7f);
		tmpno1 = _ln(tmpno3*tmpno2+1e-7f,log_table,AEC_LN_PRECISION);
		concat_buf[q] = tmpno1*a1 + a2;
	}
	//tmpno1 = logf(fabs(echo_fft_buf[0])+1e-7f);
	tmpno1 = _ln(fabs(echo_fft_buf[0])+1e-7f,log_table,AEC_LN_PRECISION);
	concat_buf[TARGET_DIM] = tmpno1*a1 + a2;
	for(i=2,j=1;i<TARGET_DIM*2;i+=2,j++)
	{
		tmpno1 = echo_fft_buf[i]*echo_fft_buf[i];
		tmpno2 = echo_fft_buf[i+1]*echo_fft_buf[i+1] + tmpno1;
		tmpno3 = _reciprocal_sqrt_hp(tmpno2+1e-8f);
		//tmpno1 = logf(sqrt(tmpno2)+1e-7f);
		tmpno1 = _ln(tmpno3*tmpno2+1e-7f,log_table,AEC_LN_PRECISION);
		concat_buf[TARGET_DIM+j] = tmpno1*a1 + a2;
	}
	//tmpno1 = logf(fabs(far_fft_buf[0])+1e-7f);
	tmpno1 = _ln(fabs(far_fft_buf[0])+1e-7f,log_table,AEC_LN_PRECISION);
	concat_buf[TARGET_DIM*2] = tmpno1*a1 + a2;
	for(m=2,n=1;m<TARGET_DIM*2;m+=2,n++)
	{
		tmpno1 = far_fft_buf[m]*far_fft_buf[m];
		tmpno2 = far_fft_buf[m+1]*far_fft_buf[m+1] + tmpno1;
		tmpno3 = _reciprocal_sqrt_hp(tmpno2+1e-8f);
		//tmpno1 = logf(sqrt(tmpno2)+1e-7f);
		tmpno1 = _ln(tmpno3*tmpno2+1e-7f,log_table,AEC_LN_PRECISION);
		concat_buf[2*TARGET_DIM+n] = tmpno1*a1 + a2;
	}

	AecPostProcess_DenseLayer_555X512_ActivationTanh(concat_buf,nn_layer_buf,exp_table,exp_precision);
	AecPostProcess_DenseLayer_512X256_ActivationTanh(nn_layer_buf,concat_buf,exp_table,exp_precision);
	AecPostProcess_DenseLayer_256X256_ActivationTanh(concat_buf,nn_layer_buf,exp_table,exp_precision);
	AecPostProcess_DenseLayer_256X185_ActivationTanh(nn_layer_buf,concat_buf,exp_table,exp_precision);

	enhanced[0] = res_fft_buf[0]*concat_buf[0];
	enhanced[1] = 0.0f;
	for(i=2,j=1;i<TARGET_DIM*2;i+=2,j++)
	{
		enhanced[i] = res_fft_buf[i]*concat_buf[j];
		enhanced[i+1] = res_fft_buf[i+1]*concat_buf[j];
	}

	memset(&enhanced[TARGET_DIM*2],0,(POST_FFT_LEN-TARGET_DIM*2)*sizeof(float));
	rdft(POST_FFT_LEN,&(fft_conf->nc),&(fft_conf->nw),-1,enhanced,fft_conf->ip,fft_conf->w);

	return 0;
}
