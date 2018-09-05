
#include "subband/analy_synth/filterbank.h"
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "aec/aec_defines.h"
#include "subband/prototype/prototype_filter.h"
#include "utility/fft/rdft_8g.h"

int DftFilterBankCreate(FilterBankControl **fb_inst) {
    FilterBankControl *fbc;
    if (fb_inst == NULL) {
        return -1;
    }

    fbc = (FilterBankControl *)malloc(sizeof(FilterBankControl));
    if (fbc == NULL) {
        return -1;
    }
    fbc->analysis_buf_head_pos = PROTOTYPE_FILTER_LEN - SUBBAND_FRAME_SHIFT;
    fbc->analysis_buf = (float *)malloc(sizeof(float) * PROTOTYPE_FILTER_LEN);
    if (fbc->analysis_buf == NULL) {
        return -1;
    }
    fbc->analysis_fft_buf = (float *)malloc(sizeof(float) * SUBBAND_NUM);
    if (fbc->analysis_fft_buf == NULL) {
        return -1;
    }
    fbc->synthesis_fft_buf = (float *)malloc(sizeof(float) * (SUBBAND_NUM + 1));
    if (fbc->synthesis_fft_buf == NULL) {
        return -1;
    }
    fbc->synthesis_state_buf =
        (float *)malloc(sizeof(float) * PROTOTYPE_FILTER_LEN);
    if (fbc->synthesis_state_buf == NULL) {
        return -1;
    }

    *fb_inst = fbc;

    return 0;
}

int DftFilterBankInit(FilterBankControl *fb_inst) {
    if (fb_inst == NULL) {
        return -1;
    }
    fb_inst->prototype_filter_len = PROTOTYPE_FILTER_LEN;
    fb_inst->oversample_rate = 2;
    fb_inst->subband_num = SUBBAND_NUM;
    fb_inst->analysis_buf_head_pos = PROTOTYPE_FILTER_LEN - SUBBAND_FRAME_SHIFT;
    memset(fb_inst->analysis_buf, 0, sizeof(float) * PROTOTYPE_FILTER_LEN);
    memset(fb_inst->analysis_fft_buf, 0, sizeof(float) * SUBBAND_NUM);
    memset(fb_inst->synthesis_fft_buf, 0, sizeof(float) * (SUBBAND_NUM + 1));
    memset(fb_inst->synthesis_state_buf, 0,
           sizeof(float) * PROTOTYPE_FILTER_LEN);

    return 0;
}

int DftFilterBankAnalysis(FilterBankControl *fb_inst, FFT_Config *fft_config,
                          float *in_buf, short sample_num) {
    short i, first_len, N1, N2, R1, R2;
    float *xk = fb_inst->analysis_fft_buf, *yk, *ptr1, *ptr2, *ptr3;
    float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;

    if (fb_inst == NULL || in_buf == NULL || sample_num != (SUBBAND_NUM >> 1) ||
        fft_config == NULL) {
        return -1;
    }
    // fb_inst->analysis_buf_head_pos points to the end of to-be-filled space
    fb_inst->analysis_buf_head_pos += SUBBAND_FRAME_SHIFT;
    if (fb_inst->analysis_buf_head_pos == PROTOTYPE_FILTER_LEN) {
        fb_inst->analysis_buf_head_pos = 0;
    }

    if (fb_inst->analysis_buf_head_pos == 0) {
        memcpy(fb_inst->analysis_buf +
                   (PROTOTYPE_FILTER_LEN - SUBBAND_FRAME_SHIFT),
               in_buf, SUBBAND_FRAME_SHIFT * sizeof(float));
    } else {
        memcpy(fb_inst->analysis_buf + fb_inst->analysis_buf_head_pos -
                   SUBBAND_FRAME_SHIFT,
               in_buf, SUBBAND_FRAME_SHIFT * sizeof(float));
    }

    yk = fb_inst->analysis_buf;
    first_len = (PROTOTYPE_FILTER_LEN - fb_inst->analysis_buf_head_pos);
    N1 = first_len / SUBBAND_FRAME_SHIFT;
    N2 = first_len / SUBBAND_NUM;
    R1 = first_len - N2 * SUBBAND_NUM;
    R2 = SUBBAND_NUM - R1;

    if (N2 == (PROTOTYPE_FILTER_LEN / SUBBAND_NUM) &&
        R1 == 0) { /*strictly SUBBAND_NUM aligned*/
        for (i = 0; i < SUBBAND_NUM; i++) {
            tmp1 = prototype_filter_coeff_44k[i] * yk[i];
            tmp2 = prototype_filter_coeff_44k[SUBBAND_NUM + i] *
                       yk[i + SUBBAND_NUM] +
                   tmp1;
            tmp3 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - i] *
                       yk[i + 2 * SUBBAND_NUM] +
                   tmp2;
            xk[i] = prototype_filter_coeff_44k[SUBBAND_NUM - 1 - i] *
                        yk[i + 3 * SUBBAND_NUM] +
                    tmp3;
        }
    } else if (N2 == 0 && R1 > 0) {
        ptr1 = &(fb_inst->analysis_buf[fb_inst->analysis_buf_head_pos]);
        ptr2 = &(fb_inst->analysis_buf[R2]);
        ptr3 = &(fb_inst->analysis_buf[0]);
        for (i = 0; i < R1; i++) {
            tmp1 = prototype_filter_coeff_44k[i] * ptr1[i];
            tmp2 = prototype_filter_coeff_44k[SUBBAND_NUM + i] * ptr2[i] + tmp1;
            tmp3 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - i] *
                       ptr2[i + SUBBAND_NUM] +
                   tmp2;
            xk[i] = prototype_filter_coeff_44k[SUBBAND_NUM - 1 - i] *
                        ptr2[i + 2 * SUBBAND_NUM] +
                    tmp3;

            tmp4 = prototype_filter_coeff_44k[i + R1] * ptr3[i];
            tmp5 = prototype_filter_coeff_44k[SUBBAND_NUM + R1 + i] *
                       ptr3[i + SUBBAND_NUM] +
                   tmp4;
            tmp6 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - R1 - i] *
                       ptr3[i + 2 * SUBBAND_NUM] +
                   tmp5;
            xk[i + R1] = prototype_filter_coeff_44k[SUBBAND_NUM - 1 - R1 - i] *
                             ptr3[i + 3 * SUBBAND_NUM] +
                         tmp6;
        }
    } else if (N2 > 0 && R1 > 0) {
        ptr1 = &(fb_inst->analysis_buf[fb_inst->analysis_buf_head_pos]);
        ptr2 = &(fb_inst->analysis_buf[R2]);
        ptr3 = &(fb_inst->analysis_buf[0]);
        switch (N2) {
            case 1: {
                for (i = 0; i < R1; i++) {
                    tmp1 = prototype_filter_coeff_44k[i] * ptr1[i];
                    tmp2 = prototype_filter_coeff_44k[SUBBAND_NUM + i] *
                               ptr1[i + SUBBAND_NUM] +
                           tmp1;
                    tmp3 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - i] *
                               ptr2[i] +
                           tmp2;
                    xk[i] = prototype_filter_coeff_44k[SUBBAND_NUM - 1 - i] *
                                ptr2[i + SUBBAND_NUM] +
                            tmp3;

                    tmp4 = prototype_filter_coeff_44k[i + R1] * ptr1[i + R1];
                    tmp5 = prototype_filter_coeff_44k[SUBBAND_NUM + R1 + i] *
                               ptr3[i] +
                           tmp4;
                    tmp6 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - R1 -
                                                      i] *
                               ptr3[i + SUBBAND_NUM] +
                           tmp5;
                    xk[i + R1] =
                        prototype_filter_coeff_44k[SUBBAND_NUM - 1 - R1 - i] *
                            ptr3[i + 2 * SUBBAND_NUM] +
                        tmp6;
                }
                break;
            }
            case 2: {
                for (i = 0; i < R1; i++) {
                    tmp1 = prototype_filter_coeff_44k[i] * ptr1[i];
                    tmp2 = prototype_filter_coeff_44k[SUBBAND_NUM + i] *
                               ptr1[i + SUBBAND_NUM] +
                           tmp1;
                    tmp3 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - i] *
                               ptr1[i + 2 * SUBBAND_NUM] +
                           tmp2;
                    xk[i] = prototype_filter_coeff_44k[SUBBAND_NUM - 1 - i] *
                                ptr2[i] +
                            tmp3;

                    tmp4 = prototype_filter_coeff_44k[i + R1] * ptr1[i + R1];
                    tmp5 = prototype_filter_coeff_44k[SUBBAND_NUM + R1 + i] *
                               ptr1[i + SUBBAND_NUM + R1] +
                           tmp4;
                    tmp6 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - R1 -
                                                      i] *
                               ptr3[i] +
                           tmp5;
                    xk[i + R1] =
                        prototype_filter_coeff_44k[SUBBAND_NUM - 1 - R1 - i] *
                            ptr3[i + SUBBAND_NUM] +
                        tmp6;
                }
                break;
            }
            case 3: {
                for (i = 0; i < R1; i++) {
                    tmp1 = prototype_filter_coeff_44k[i] * ptr1[i];
                    tmp2 = prototype_filter_coeff_44k[SUBBAND_NUM + i] *
                               ptr1[i + SUBBAND_NUM] +
                           tmp1;
                    tmp3 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - i] *
                               ptr1[i + 2 * SUBBAND_NUM] +
                           tmp2;
                    xk[i] = prototype_filter_coeff_44k[SUBBAND_NUM - 1 - i] *
                                ptr1[i + 3 * SUBBAND_NUM] +
                            tmp3;

                    tmp4 = prototype_filter_coeff_44k[i + R1] * ptr1[i + R1];
                    tmp5 = prototype_filter_coeff_44k[SUBBAND_NUM + R1 + i] *
                               ptr1[i + SUBBAND_NUM + R1] +
                           tmp4;
                    tmp6 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - R1 -
                                                      i] *
                               ptr1[i + 2 * SUBBAND_NUM + R1] +
                           tmp5;
                    xk[i + R1] =
                        prototype_filter_coeff_44k[SUBBAND_NUM - 1 - R1 - i] *
                            ptr3[i] +
                        tmp6;
                }
                break;
            }
            default:
                break;
        }
    } else if (N2 > 0 && R1 == 0) {
        ptr1 = &(fb_inst->analysis_buf[fb_inst->analysis_buf_head_pos]);
        ptr2 = &(fb_inst->analysis_buf[0]);
        switch (N2) {
            case 1: {
                for (i = 0; i < SUBBAND_NUM; i++) {
                    tmp1 = prototype_filter_coeff_44k[i] * ptr1[i];
                    tmp2 =
                        prototype_filter_coeff_44k[SUBBAND_NUM + i] * ptr2[i] +
                        tmp1;
                    tmp3 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - i] *
                               ptr2[i + SUBBAND_NUM] +
                           tmp2;
                    xk[i] = prototype_filter_coeff_44k[SUBBAND_NUM - 1 - i] *
                                ptr2[i + 2 * SUBBAND_NUM] +
                            tmp3;
                }
                break;
            }
            case 2: {
                for (i = 0; i < SUBBAND_NUM; i++) {
                    tmp1 = prototype_filter_coeff_44k[i] * ptr1[i];
                    tmp2 = prototype_filter_coeff_44k[SUBBAND_NUM + i] *
                               ptr1[SUBBAND_NUM + i] +
                           tmp1;
                    tmp3 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - i] *
                               ptr2[i] +
                           tmp2;
                    xk[i] = prototype_filter_coeff_44k[SUBBAND_NUM - 1 - i] *
                                ptr2[i + SUBBAND_NUM] +
                            tmp3;
                }
                break;
            }
            case 3: {
                for (i = 0; i < SUBBAND_NUM; i++) {
                    tmp1 = prototype_filter_coeff_44k[i] * ptr1[i];
                    tmp2 = prototype_filter_coeff_44k[SUBBAND_NUM + i] *
                               ptr1[SUBBAND_NUM + i] +
                           tmp1;
                    tmp3 = prototype_filter_coeff_44k[2 * SUBBAND_NUM - 1 - i] *
                               ptr1[2 * SUBBAND_NUM + i] +
                           tmp2;
                    xk[i] = prototype_filter_coeff_44k[SUBBAND_NUM - 1 - i] *
                                ptr2[i] +
                            tmp3;
                }
                break;
            }
            default:
                break;
        }
    }

    rdft(SUBBAND_NUM, &(fft_config->nc), &(fft_config->nw), 1, xk,
         fft_config->ip, fft_config->w);

    return 0;
}

int DftFilterBankSynthesis(FilterBankControl *fb_inst, FFT_Config *fft_config,
                           float *out_buf, short *sample_num) {
    short D = SUBBAND_NUM / OVERSAMPLE_RATE, i;
    float *buf1, *buf2, *buf3, tmp1, tmp2, gain;
    float *flt_seg1, *flt_seg2, *flt_seg3, *flt_seg4, *flt_seg5, *flt_seg6,
        *flt_seg7;

    if (fb_inst == NULL || out_buf == NULL || fft_config == NULL) {
        return -1;
    }
    /*here rdft is used instead of cdft, so the sequence order and gain need to
     * adjust accordingly*/
    rdft(SUBBAND_NUM, &(fft_config->nc), &(fft_config->nw), -1,
         fb_inst->synthesis_fft_buf, fft_config->ip, fft_config->w);

    buf1 = &(fb_inst->synthesis_fft_buf[D]);
    buf2 = &(fb_inst->synthesis_state_buf[PROTOTYPE_FILTER_LEN - D]);
    buf3 = &(prototype_filter_coeff_44k[D - 1]);
    gain = (PROTOTYPE_FILTER_LEN == 2048) ? 256.0f : 128.0f;
    for (i = 0; i < D; i++) {
        tmp1 = *buf1-- * 2.0f;
        tmp2 = tmp1 * *buf3-- + *buf2++;
        out_buf[D - 1 - i] = tmp2 * gain;
    }
    *sample_num = D;
    /*matlab DFT filter-bank synthesis, in = (repmat(A(:,k),mlp,1))';v = v +
     * h.*in;*/
    /*divide filter into 7 segments*/
    flt_seg1 = &(prototype_filter_coeff_44k[0]);
    flt_seg2 = &(prototype_filter_coeff_44k[PROTOTYPE_FILTER_LEN >> 3]);
    flt_seg3 = &(prototype_filter_coeff_44k[PROTOTYPE_FILTER_LEN >> 2]);
    flt_seg4 = &(prototype_filter_coeff_44k[(PROTOTYPE_FILTER_LEN >> 2) +
                                            (PROTOTYPE_FILTER_LEN >> 3)]);
    flt_seg5 = &(prototype_filter_coeff_44k[(PROTOTYPE_FILTER_LEN >> 1) - 1]);
    flt_seg6 = &(prototype_filter_coeff_44k[(PROTOTYPE_FILTER_LEN >> 2) +
                                            (PROTOTYPE_FILTER_LEN >> 3) - 1]);
    flt_seg7 = &(prototype_filter_coeff_44k[(PROTOTYPE_FILTER_LEN >> 2) - 1]);

    /*divide synthesis_fft_buf into 2 segments*/
    fb_inst->synthesis_fft_buf[SUBBAND_NUM] = fb_inst->synthesis_fft_buf[0];
    buf1 = &(fb_inst->synthesis_fft_buf[SUBBAND_NUM]);
    buf2 = &(fb_inst->synthesis_fft_buf[SUBBAND_NUM / 2]);

    // buf3 = &(fb_inst->synthesis_state_buf[2*PROTOTYPE_FILTER_LEN-D])
    buf3 = &(fb_inst->synthesis_state_buf[0]);
    for (i = 0; i < D; i++) {
        tmp1 = *buf1-- * 2.0f;
        tmp2 = *buf2-- * 2.0f;
        buf3[PROTOTYPE_FILTER_LEN - 1 * D + i] =
            buf3[PROTOTYPE_FILTER_LEN - 2 * D + i] + *flt_seg7-- * tmp1;
        buf3[PROTOTYPE_FILTER_LEN - 2 * D + i] =
            buf3[PROTOTYPE_FILTER_LEN - 3 * D + i] + *flt_seg6-- * tmp2;
        buf3[PROTOTYPE_FILTER_LEN - 3 * D + i] =
            buf3[PROTOTYPE_FILTER_LEN - 4 * D + i] + *flt_seg5-- * tmp1;
        buf3[PROTOTYPE_FILTER_LEN - 4 * D + i] =
            buf3[PROTOTYPE_FILTER_LEN - 5 * D + i] + *flt_seg4++ * tmp2;
        buf3[PROTOTYPE_FILTER_LEN - 5 * D + i] =
            buf3[PROTOTYPE_FILTER_LEN - 6 * D + i] + *flt_seg3++ * tmp1;
        buf3[PROTOTYPE_FILTER_LEN - 6 * D + i] =
            buf3[PROTOTYPE_FILTER_LEN - 7 * D + i] + *flt_seg2++ * tmp2;
        buf3[PROTOTYPE_FILTER_LEN - 7 * D + i] =
            buf3[PROTOTYPE_FILTER_LEN - 8 * D + i] + *flt_seg1++ * tmp1;
        buf3[i] = 0.0f;
    }

    return 0;
}

int DftFilterBankReset(FilterBankControl *fb_inst) {
    if (fb_inst == NULL) {
        return -1;
    }

    fb_inst->prototype_filter_len = PROTOTYPE_FILTER_LEN;
    fb_inst->oversample_rate = 2;
    fb_inst->subband_num = SUBBAND_NUM;
    fb_inst->analysis_buf_head_pos = PROTOTYPE_FILTER_LEN - SUBBAND_FRAME_SHIFT;
    memset(fb_inst->analysis_buf, 0, sizeof(float) * PROTOTYPE_FILTER_LEN);
    memset(fb_inst->analysis_fft_buf, 0, sizeof(float) * SUBBAND_NUM);
    memset(fb_inst->synthesis_fft_buf, 0, sizeof(float) * (SUBBAND_NUM + 1));
    memset(fb_inst->synthesis_state_buf, 0,
           sizeof(float) * PROTOTYPE_FILTER_LEN);

    return 0;
}

int DftFilterBankFree(FilterBankControl *fb_inst) {
    if (fb_inst == NULL) {
        return -1;
    }

    if (fb_inst->analysis_buf != NULL) {
        free(fb_inst->analysis_buf);
    }

    if (fb_inst->analysis_fft_buf != NULL) {
        free(fb_inst->analysis_fft_buf);
    }

    if (fb_inst->synthesis_fft_buf != NULL) {
        free(fb_inst->synthesis_fft_buf);
    }

    if (fb_inst->synthesis_state_buf != NULL) {
        free(fb_inst->synthesis_state_buf);
    }

    free(fb_inst);

    return 0;
}
