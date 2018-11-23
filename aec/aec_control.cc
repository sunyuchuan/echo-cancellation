#ifdef __cplusplus
extern "C" {
#endif

#include "aec/aec_control.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aec/adp_filter_coeff_factor.h"
#include "aec/aec_core.h"
#include "aec/aec_defines.h"
#include "aec/delay_estimator/delay_estimator_internal.h"
#include "aec/delay_estimator/delay_estimator_wrapper.h"
#include "subband/analy_synth/filterbank.h"
#include "subband/prototype/prototype_defines.h"
#include "utility/fft/fft_wrapper.h"
#include "utility/fft/rdft_8g_init.h"
#include "utility/log/log.h"

#define IFFT_GAIN (32767.0f * 2.0f / 1024.0f)

namespace xmly_audio_recorder_android {
AecControl::AecControl()
    : far_buf_write_pos(0),
      far_buf_read_pos(0),
      known_delay(0),
      last_known_delay(0),
      far_frame_buf(NULL),
      near_frame_buf(NULL),
      far_end_buf(NULL),
      last_far_end_frame(NULL),
      delay_estimator_farend(NULL),
      delay_estimator(NULL),
      current_delay(-2),
      far_history(NULL),
      far_sepctrum_history(NULL),
      far_history_pos(0),
      realFFT(NULL),
      fb_ctrl_far(NULL),
      fb_ctrl_near(NULL),
      near_input_buf(NULL),
      far_input_buf(NULL),
      adp_filter_coeff(NULL),
      far_fb_buf(NULL),
      near_fb_buf(NULL),
      echo_buf(NULL),
      error_buf(NULL),
      post_res_frame(NULL),
#if AEC_POST_PROCESSING_ON
      post_fft_conf(NULL),
      coh_res(NULL),
      coh_echo(NULL),
      coh_rd(NULL),
      adpt_coh_echo(NULL),
      prev_enchanced_sqrd(NULL),
      enhanced(NULL),
      nonlinearity(1.0f),
      first_frame(1),
      post_echo_frame(NULL),
      vlp_buf(NULL),
#endif
      Rss(NULL),
      Rdd(NULL),
      Ree(NULL) {
}

int AecControl::AudioProcessing_AEC_Create() {
    far_frame_buf = (short *)malloc((FRAME_LEN + PART_LEN) * sizeof(short));
    if (far_frame_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    near_frame_buf = (short *)malloc((FRAME_LEN + PART_LEN) * sizeof(short));
    if (near_frame_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    delay_estimator_farend =
        DelayEstimator_CreateDelayEstimatorFarend(PART_LEN1, MAX_DELAY);
    if (delay_estimator_farend == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    delay_estimator =
        DelayEstimator_CreateDelayEstimator(delay_estimator_farend, 0);
    if (delay_estimator == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    far_history = (float *)malloc(sizeof(float) * PART_LEN1 * MAX_DELAY);
    if (far_history == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    far_sepctrum_history =
        (float *)malloc(sizeof(float) * SUBBAND_NUM * MAX_DELAY);
    if (far_sepctrum_history == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    /*initialize fft*/
    if (rdft_create(&realFFT, SUBBAND_NUM) < 0) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    /*initialize filter-bank*/
    if (DftFilterBankCreate(&fb_ctrl_far) < 0) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    if (DftFilterBankCreate(&fb_ctrl_near) < 0) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    near_input_buf = (float *)malloc(sizeof(float) * PROTOTYPE_FILTER_LEN);
    if (near_input_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    far_input_buf = (float *)malloc(sizeof(float) * PROTOTYPE_FILTER_LEN);
    if (far_input_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    /*oversample rate is OVERSAMPLE_RATE and 1 complex number consists of 2
     * float numbers so ADPF_LEN*OVERSAMPLE_RATE*2*/
    adp_filter_coeff =
        (float *)malloc(sizeof(float) * (ADPF_LEN * OVERSAMPLE_RATE) * 2);
    if (adp_filter_coeff == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }
    memcpy(adp_filter_coeff, adp_filter_coeff_factor,
           sizeof(float) * ADP_FILTER_COEFF_FACTOR_LEN);

    far_fb_buf =
        (float *)malloc(sizeof(float) * (ADPF_LEN * OVERSAMPLE_RATE) * 2);
    if (far_fb_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    near_fb_buf =
        (float *)malloc(sizeof(float) * (ADPF_LEN * OVERSAMPLE_RATE) * 2);
    if (near_fb_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    echo_buf = (float *)malloc(sizeof(float) * 2 * (SUBBAND_NUM / 2));
    if (echo_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    error_buf = (float *)malloc(sizeof(float) * 2 * (SUBBAND_NUM / 2));
    if (error_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    far_end_buf = RingBuffer_CreateBuffer(
        SUBBAND_FRAME_SHIFT * FAR_END_BUF_FRAME_NUM, sizeof(short));
    if (far_end_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    last_far_end_frame = (short *)malloc(sizeof(short) * FRAME_LEN);
    if (last_far_end_frame == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    post_res_frame = (float *)malloc(sizeof(float) * (POST_FFT_LEN));
    if (post_res_frame == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }
#if AEC_POST_PROCESSING_ON
    if (rdft_create(&post_fft_conf, POST_FFT_LEN) < 0) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    coh_res = (float *)malloc(sizeof(float) * (POST_FFT_LEN / 2 + 1));
    if (coh_res == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    coh_echo = (float *)malloc(sizeof(float) * (POST_FFT_LEN / 2 + 1));
    if (coh_echo == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    coh_rd = (float *)malloc(sizeof(float) * (POST_FFT_LEN));
    if (coh_echo == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    adpt_coh_echo = (float *)malloc(sizeof(float) * (POST_FFT_LEN / 2 + 1));
    if (adpt_coh_echo == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    prev_enchanced_sqrd =
        (float *)malloc(sizeof(float) * (POST_FFT_LEN / 2 + 1));
    if (prev_enchanced_sqrd == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    enhanced = (float *)malloc(sizeof(float) * (POST_FFT_LEN));
    if (enhanced == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    post_echo_frame = (float *)malloc(sizeof(float) * (POST_FFT_LEN));
    if (post_echo_frame == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    vlp_buf = (float *)malloc(sizeof(float) * (POST_FFT_LEN / 2));
    if (vlp_buf == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }
#endif

    Rss = (float *)malloc(sizeof(float) * (SUBBAND_NUM / 2 + 1));
    if (Rss == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    Rdd = (float *)malloc(sizeof(float) * (SUBBAND_NUM / 2 + 1));
    if (Rdd == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    Ree = (float *)malloc(sizeof(float) * (SUBBAND_NUM / 2 + 1));
    if (Ree == NULL) {
        this->AudioProcessing_AEC_Release();
        return -1;
    }

    return 0;
}

int AecControl::AudioProcessing_AEC_Init(float amp_perc, float min_perc) {
    far_buf_write_pos = 0;
    far_buf_read_pos = 0;
    known_delay = 0;
    last_known_delay = 0;

    SetNonlinearGain(amp_perc, &nonlinearity, min_perc);

    RingBuffer_InitBuffer(far_end_buf);

    memset(far_frame_buf, 0, (FRAME_LEN + PART_LEN) * sizeof(short));
    memset(near_frame_buf, 0, (FRAME_LEN + PART_LEN) * sizeof(short));
    memset(last_far_end_frame, 0, sizeof(short) * FRAME_LEN);
    memcpy(adp_filter_coeff, adp_filter_coeff_factor,
           sizeof(float) * ADP_FILTER_COEFF_FACTOR_LEN);
    memset(far_fb_buf, 0, sizeof(float) * (ADPF_LEN * OVERSAMPLE_RATE) * 2);
    memset(near_fb_buf, 0, sizeof(float) * (ADPF_LEN * OVERSAMPLE_RATE) * 2);
    memset(error_buf, 0, sizeof(float) * 2 * (SUBBAND_NUM / 2));
    memset(echo_buf, 0, sizeof(float) * 2 * (SUBBAND_NUM / 2));
    memset(Rss, 0, sizeof(float) * (SUBBAND_NUM / 2 + 1));
    memset(Rdd, 0, sizeof(float) * (SUBBAND_NUM / 2 + 1));
    memset(Ree, 0, sizeof(float) * (SUBBAND_NUM / 2 + 1));
    memset(post_res_frame, 0, sizeof(float) * (POST_FFT_LEN));
#if AEC_POST_PROCESSING_ON
    memset(coh_res, 0, sizeof(float) * (POST_FFT_LEN / 2 + 1));
    memset(coh_echo, 0, sizeof(float) * (POST_FFT_LEN / 2 + 1));
    memset(coh_rd, 0, sizeof(float) * (POST_FFT_LEN));
    memset(adpt_coh_echo, 0, sizeof(float) * (POST_FFT_LEN / 2 + 1));
    memset(prev_enchanced_sqrd, 0, sizeof(float) * (POST_FFT_LEN / 2 + 1));
    memset(enhanced, 0, sizeof(float) * (POST_FFT_LEN));
    memset(post_echo_frame, 0, sizeof(float) * (POST_FFT_LEN));
    memset(vlp_buf, 0, sizeof(float) * (POST_FFT_LEN / 2));
    rdft_init(post_fft_conf->ip, post_fft_conf->w, POST_FFT_LEN,
              &post_fft_conf->nw, &post_fft_conf->nc);
#endif

    rdft_init(realFFT->ip, realFFT->w, (PART_LEN << 1), &realFFT->nw,
              &realFFT->nc);

    DftFilterBankInit(fb_ctrl_far);
    DftFilterBankInit(fb_ctrl_near);

    if (DelayEstimator_InitDelayEstimatorFarend(delay_estimator_farend) != 0) {
        return -1;
    }
    if (DelayEstimator_InitDelayEstimator(delay_estimator) != 0) {
        return -1;
    }
    memset(far_history, 0, sizeof(float) * PART_LEN1 * MAX_DELAY);
    memset(far_sepctrum_history, 0, sizeof(float) * SUBBAND_NUM * MAX_DELAY);
    far_history_pos = MAX_DELAY;

    return 0;
}

int AecControl::AudioProcessing_AEC_FillFarBuf(char *ref, short sample_size,
                                               bool playout_switch,
                                               bool mic_switch) {
    if (ref == NULL || sample_size < 0) {
        // LOGI("input arguments error.\n");
        return -1;
    }

    if ((sample_size % SUBBAND_FRAME_SHIFT) != 0) {
        // LOGI("sample_size invalid.\n");
        return -1;
    }

    if (playout_switch == true) {
        RingBuffer_WriteBuffer(far_end_buf, ref, (size_t)(sample_size >> 1));
    }

    if (playout_switch == true && mic_switch == false) {
        short farend[FRAME_LEN], i, k;
        const short *farend_ptr = NULL;
        float freq_dmn_far_inv_abs[SUBBAND_NUM / 2 + 1],
            freq_dmn_far_abs[SUBBAND_NUM / 2 + 1];
        short num_of_filled_buffers =
            (short)RingBuffer_available_read(far_end_buf) / SUBBAND_FRAME_SHIFT;

        if (num_of_filled_buffers >= 16) {
            for (i = 0; i < 16; i++) {
                RingBuffer_ReadBuffer(far_end_buf, (void **)&farend_ptr, farend,
                                      SUBBAND_FRAME_SHIFT);
                memcpy(last_far_end_frame, farend_ptr,
                       SUBBAND_FRAME_SHIFT * sizeof(short));

                for (k = 0; k < SUBBAND_FRAME_SHIFT; k++) {
                    far_input_buf[k] = ((float)farend_ptr[k]) * 3.05185094e-5f;
                    // near_input_buf[k] = ((float)pri[k])*3.05185094e-5f;
                }

                DftFilterBankAnalysis(fb_ctrl_far, realFFT, far_input_buf,
                                      SUBBAND_FRAME_SHIFT);
                CalcAbsValue(fb_ctrl_far, freq_dmn_far_inv_abs,
                             freq_dmn_far_abs);
                UpdateFarHistory(&far_history_pos, far_history,
                                 freq_dmn_far_abs, far_sepctrum_history,
                                 fb_ctrl_far->analysis_fft_buf);
                if (DelayEstimator_AddFarSpectrumFloat(delay_estimator_farend,
                                                       freq_dmn_far_abs,
                                                       PART_LEN1) < 0) {
                    return -1;
                }
            }
        }
    }
    // LOGI("far_end_buf->read_pos = %lu, far_end_buf->write_pos = %lu.\n",
    //      far_end_buf->read_pos, far_end_buf->write_pos);

    return 0;
}

int AecControl::AudioProcessing_AEC_Process(
    char *pri, short sample_size, char *out_buf, unsigned int *output_size,
    char *out_buf1, unsigned int *output_size1, bool playout_switch,
    bool mic_switch, bool amp_pwr_changed, float amp_level, float min_level) {
    /*if(pri==NULL || sample_size<0)
    {
            return -1;
    }*/
    short i, k, q;
    short farend[FRAME_LEN];
    const short *farend_ptr = NULL;
    float freq_dmn_far_inv_abs[SUBBAND_NUM / 2 + 1],
        freq_dmn_near_inv_abs[SUBBAND_NUM / 2 + 1],
        freq_dmn_far_abs[SUBBAND_NUM / 2 + 1],
        freq_dmn_near_abs[SUBBAND_NUM / 2 + 1], output[SUBBAND_NUM / 2],
        output1[SUBBAND_NUM / 2];

    short nFrames = (sample_size >> 1) / SUBBAND_FRAME_SHIFT, output_sample_num,
          output_sample_num1, seg_len = 0, data_len_record1 = POST_FFT_LEN / 2,
          data_len_record2 = POST_FFT_LEN / 2, *out = (short *)out_buf,
          *near = (short *)pri;
    float *far_spectrum, tmpno1, tmpno2, tmpno3, tmpno4;
    int delay;

    if (((sample_size >> 1) - nFrames * SUBBAND_FRAME_SHIFT) !=
        0) { /*num_of_samples should be exactly integer multiplies of
                SUBBAND_FRAME_SHIFT*/
        return -1;
    }

    *output_size1 = 0;
    *output_size = 0;
    if (amp_pwr_changed) {
        SetNonlinearGain(amp_level, &nonlinearity, min_level);
    }

    for (i = 0; i < nFrames; i++) {
        if (playout_switch == true) {
            /*check if there is enough data in far end buffer*/
            short num_of_filled_buffers =
                (short)RingBuffer_available_read(far_end_buf) /
                SUBBAND_FRAME_SHIFT;

            if (num_of_filled_buffers > 0) {
                RingBuffer_ReadBuffer(far_end_buf, (void **)&farend_ptr, farend,
                                      SUBBAND_FRAME_SHIFT);
                memcpy(last_far_end_frame, farend_ptr,
                       SUBBAND_FRAME_SHIFT * sizeof(short));
            } else { /*use last buffer instead if far end can't provide any
                        data*/
                memcpy(farend, last_far_end_frame,
                       SUBBAND_FRAME_SHIFT * sizeof(short));
                farend_ptr = farend;
            }
            for (k = 0; k < SUBBAND_FRAME_SHIFT; k++) {
                far_input_buf[k] = ((float)farend_ptr[k]) * 3.05185094e-5f;
                // near_input_buf[k] = ((float)pri[k])*3.05185094e-5f;
            }

            DftFilterBankAnalysis(fb_ctrl_far, realFFT, far_input_buf,
                                  SUBBAND_FRAME_SHIFT);
            CalcAbsValue(fb_ctrl_far, freq_dmn_far_inv_abs, freq_dmn_far_abs);
            UpdateFarHistory(&far_history_pos, far_history, freq_dmn_far_abs,
                             far_sepctrum_history,
                             fb_ctrl_far->analysis_fft_buf);
        }

        if (mic_switch == true) {
            for (k = 0; k < SUBBAND_FRAME_SHIFT; k++) {
                // far_input_buf[k] = ((float)farend_ptr[k])*3.05185094e-5f;
                near_input_buf[k] =
                    ((float)near[i * SUBBAND_FRAME_SHIFT + k]) * 3.05185094e-5f;
            }
            DftFilterBankAnalysis(fb_ctrl_near, realFFT, near_input_buf,
                                  SUBBAND_FRAME_SHIFT);
            CalcAbsValue(fb_ctrl_near, freq_dmn_near_inv_abs,
                         freq_dmn_near_abs);
        }

        if (mic_switch == true && playout_switch == true) {
            if (DelayEstimator_AddFarSpectrumFloat(
                    delay_estimator_farend, freq_dmn_far_abs, PART_LEN1) < 0) {
                return -1;
            }

            delay = DelayEstimator_DelayEstimatorProcessFloat(
                delay_estimator, freq_dmn_near_abs, PART_LEN1);
            if (delay == -2) {
                delay = 0;
            } else if (delay == -1) {
                return -1;
            }
            // LOGI("Estimated delay is %d\n", delay);

            far_spectrum =
                AlignedFarend(far_history_pos, far_sepctrum_history, delay);

            // 5.linear echo cancellation
            AecDeecho(far_spectrum, fb_ctrl_near->analysis_fft_buf, far_fb_buf,
                      near_fb_buf, adp_filter_coeff,
                      fb_ctrl_far->synthesis_fft_buf,
                      fb_ctrl_near->synthesis_fft_buf, freq_dmn_near_abs, Rss,
                      Rdd, Ree);
        } else if (mic_switch == true && playout_switch == false) {
            memcpy(fb_ctrl_near->synthesis_fft_buf,
                   fb_ctrl_near->analysis_fft_buf, sizeof(float) * SUBBAND_NUM);
        }

        if (mic_switch == true) {
            // memcpy(fb_ctrl_near->synthesis_fft_buf,fb_ctrl_near->analysis_fft_buf,sizeof(float)*SUBBAND_NUM);
            DftFilterBankSynthesis(fb_ctrl_near, realFFT, output,
                                   &output_sample_num);
            *output_size += output_sample_num << 1;

            memcpy(post_res_frame + data_len_record1, output,
                   output_sample_num * sizeof(float));
            data_len_record1 += output_sample_num;
#if AEC_POST_PROCESSING_ON
            if (playout_switch == true) {
                DftFilterBankSynthesis(fb_ctrl_far, realFFT, output1,
                                       &output_sample_num1);
                *output_size1 += output_sample_num1 << 1;
                memcpy(post_echo_frame + data_len_record2, output1,
                       output_sample_num * sizeof(float));
                data_len_record2 += output_sample_num;
            }

            if (data_len_record1 == POST_FFT_LEN && playout_switch == true) {
                AecResidualEchoCancellation(
                    post_res_frame, post_echo_frame, post_fft_conf, coh_res,
                    coh_echo, coh_rd, adpt_coh_echo, prev_enchanced_sqrd,
                    enhanced, &first_frame, nonlinearity);
                memmove(post_res_frame, post_res_frame + (POST_FFT_LEN >> 1),
                        sizeof(float) * (POST_FFT_LEN >> 1));
                memmove(post_echo_frame, post_echo_frame + (POST_FFT_LEN >> 1),
                        sizeof(float) * (POST_FFT_LEN >> 1));
                data_len_record1 = POST_FFT_LEN >> 1;
                data_len_record2 = POST_FFT_LEN >> 1;
                for (q = 0; q < (POST_FFT_LEN >> 1); q += 2) {
                    tmpno1 = vlp_buf[q] + enhanced[q];
                    tmpno2 = vlp_buf[q + 1] + enhanced[q + 1];
                    tmpno3 = fmax(tmpno1 * IFFT_GAIN, -32768.0f);
                    tmpno4 = fmax(tmpno2 * IFFT_GAIN, -32768.0f);
                    tmpno1 = fmin(tmpno3, 32767.0f);
                    tmpno2 = fmin(tmpno4, 32767.0f);
                    out[seg_len + q] = (short)(tmpno1);
                    out[seg_len + q + 1] = (short)(tmpno2);
                    vlp_buf[q] = enhanced[POST_FFT_LEN / 2 + q];
                    vlp_buf[q + 1] = enhanced[POST_FFT_LEN / 2 + q + 1];
                }
                seg_len += (POST_FFT_LEN >> 1);
            } else
#endif
                if (data_len_record1 == POST_FFT_LEN) {
                for (q = 0; q < (POST_FFT_LEN >> 1); q += 2) {
                    tmpno1 = fmax(post_res_frame[q] * 32768.0f, -32768.0f);
                    tmpno2 = fmax(post_res_frame[q + 1] * 32768.0f, -32768.0f);
                    tmpno3 = fmin(tmpno1, 32767.0f);
                    tmpno4 = fmin(tmpno2, 32767.0f);
                    out[seg_len + q] = (short)(tmpno3);
                    out[seg_len + q + 1] = (short)(tmpno4);
                }
                memmove(post_res_frame, post_res_frame + (POST_FFT_LEN >> 1),
                        sizeof(float) * (POST_FFT_LEN >> 1));
                data_len_record1 = POST_FFT_LEN >> 1;
                seg_len += (POST_FFT_LEN >> 1);
            }
        }
    }

    return 0;
}

int AecControl::AudioProcessing_AEC_ClearFarFrameBuf() {
    memset(far_fb_buf, 0, sizeof(float) * (ADPF_LEN * OVERSAMPLE_RATE) * 2);

    return 0;
}

int AecControl::AudioProcessing_AEC_ResetNearState() {
    DftFilterBankReset(fb_ctrl_near);
#if AEC_POST_PROCESSING_ON
    memset(post_echo_frame, 0, sizeof(float) * (POST_FFT_LEN));
    memset(post_res_frame, 0, sizeof(float) * (POST_FFT_LEN));
    memset(prev_enchanced_sqrd, 0, sizeof(float) * (POST_FFT_LEN / 2 + 1));
    memset(vlp_buf, 0, sizeof(float) * (POST_FFT_LEN / 2));
    first_frame = 1;
#endif

    return 0;
}

int AecControl::AudioProcessing_AEC_ResetFarState() {
    short i;
    DelayEstimatorFarend *delay_far =
        (DelayEstimatorFarend *)delay_estimator_farend;
    memset(delay_far->mean_far_spectrum, 0,
           sizeof(SpectrumType) * delay_far->spectrum_size);
    delay_far->far_spectrum_initialized = 0;
    memset(delay_far->binary_farend->binary_far_history, 0,
           sizeof(unsigned int) * delay_far->binary_farend->history_size);
    memset(delay_far->binary_farend->far_bit_counts, 0,
           sizeof(int) * delay_far->binary_farend->history_size);

    DelayEstimator *delay_near = (DelayEstimator *)delay_estimator;
    memset(delay_near->mean_near_spectrum, 0,
           sizeof(SpectrumType) * delay_near->spectrum_size);
    delay_near->near_spectrum_initialized = 0;
    for (i = 0; i < delay_near->binary_handle->farend->history_size; ++i) {
        delay_near->binary_handle->mean_bit_counts[i] = (20 << 9);  // 20 in Q9.
    }
    delay_near->binary_handle->minimum_probability = (32 << 9);     // 32 in Q9.
    delay_near->binary_handle->last_delay_probability = (32 << 9);  // 32 in Q9.

    // Default return value if it is unable to estimate. -1 is used for errors.
    delay_near->binary_handle->last_delay = -2;

    far_end_buf->read_pos = 0;
    far_end_buf->write_pos = 0;
    far_end_buf->rw_wrap = SAME_WRAP;
    memset(far_end_buf->data, 0,
           SUBBAND_FRAME_SHIFT * FAR_END_BUF_FRAME_NUM * sizeof(short));

    far_history_pos = MAX_DELAY;
    memset(far_history, 0, sizeof(float) * PART_LEN1 * MAX_DELAY);
    memset(far_sepctrum_history, 0, sizeof(float) * SUBBAND_NUM * MAX_DELAY);

    DftFilterBankReset(fb_ctrl_far);

    return 0;
}

int AecControl::AudioProcessing_AEC_Release() {
    if (far_frame_buf != NULL) {
        free(far_frame_buf);
    }

    if (near_frame_buf != NULL) {
        free(near_frame_buf);
    }

    RingBuffer_FreeBuffer(far_end_buf);

    if (last_far_end_frame != NULL) {
        free(last_far_end_frame);
    }

    DelayEstimator_FreeDelayEstimator(delay_estimator);
    DelayEstimator_FreeDelayEstimatorFarend(delay_estimator_farend);

    if (far_history != NULL) {
        free(far_history);
    }

    if (far_sepctrum_history != NULL) {
        free(far_sepctrum_history);
    }

    if (rdft_free(realFFT) < 0) {
        return -1;
    }

    if (DftFilterBankFree(fb_ctrl_far) < 0) {
        return -1;
    }

    if (DftFilterBankFree(fb_ctrl_near) < 0) {
        return -1;
    }

    if (near_input_buf != NULL) {
        free(near_input_buf);
    }

    if (far_input_buf != NULL) {
        free(far_input_buf);
    }

    if (adp_filter_coeff != NULL) {
        free(adp_filter_coeff);
        adp_filter_coeff = NULL;
    }

    if (far_fb_buf != NULL) {
        free(far_fb_buf);
    }

    if (near_fb_buf != NULL) {
        free(near_fb_buf);
    }

    if (echo_buf != NULL) {
        free(echo_buf);
    }

    if (error_buf != NULL) {
        free(error_buf);
    }

    if (Rss != NULL) {
        free(Rss);
    }

    if (Rdd != NULL) {
        free(Rdd);
    }

    if (Ree != NULL) {
        free(Ree);
    }

    if (post_res_frame != NULL) {
        free(post_res_frame);
    }
#if AEC_POST_PROCESSING_ON
    if (rdft_free(post_fft_conf) < 0) {
        return -1;
    }

    if (coh_res != NULL) {
        free(coh_res);
    }

    if (coh_echo != NULL) {
        free(coh_echo);
    }

    if (coh_rd != NULL) {
        free(coh_rd);
    }

    if (adpt_coh_echo != NULL) {
        free(adpt_coh_echo);
    }

    if (prev_enchanced_sqrd != NULL) {
        free(prev_enchanced_sqrd);
    }

    if (enhanced != NULL) {
        free(enhanced);
    }

    if (post_echo_frame != NULL) {
        free(post_echo_frame);
    }

    if (vlp_buf != NULL) {
        free(vlp_buf);
    }

#endif

    delete this;
    return 0;
}
}  // namespace xmly_audio_recorder_android

#ifdef __cplusplus
}
#endif
