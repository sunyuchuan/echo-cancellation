#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aec/aec_control.h"
#include "aec/aec_core.h"
#include "aec/aec_defines.h"
#include "aec/delay_estimator/delay_estimator.h"
#include "subband/analy_synth/filterbank.h"

using namespace xmly_audio_recorder_android;
int main() {
    AecControl *aec_inst = new AecControl;
    FILE *fid_res, *fid_echo, *fid_out;
    short *res, *echo, out[2048], out1[2048];
    short block_size = 2048, frame_len = 1024, frame_shift = 512, i, q;
    float *res_float, *echo_float, gain = 32767.0f * 2.0f / 1024.0f, tmpno1,
                                   tmpno2;
    unsigned int size, size1, count = 0;
    bool mic, playout;
    short reset_far_state = 0, reset_near_state = 0;

    fid_res = fopen("error_pri11-aligned.pcm", "rb+");
    fid_echo = fopen("echo_pri11-aligned.pcm", "rb+");
    fid_out = fopen("out.pcm", "wb+");

    res = (short *)malloc(frame_len);
    echo = (short *)malloc(frame_len);
    res_float = (float *)malloc(sizeof(float) * frame_shift);
    echo_float = (float *)malloc(sizeof(float) * frame_shift);

    if (aec_inst == NULL) {
        return -1;
    }

    if (aec_inst->AudioProcessing_AEC_Create() < 0) {
        return -1;
    }
    if (aec_inst->AudioProcessing_AEC_Init() < 0) {
        return -1;
    }
    if (OpenDelayRecordFile() < 0) {
        return -1;
    }
    mic = true;
    playout = true;

    while (fread(res, sizeof(short), frame_shift, fid_res) == frame_shift &&
           fread(echo, sizeof(short), frame_shift, fid_echo) == frame_shift) {
        count++;
        for (i = 0; i < frame_shift; i++) {
            res_float[i] = ((float)res[i]) * 3.05185094e-5f;
            echo_float[i] = ((float)echo[i]) * 3.05185094e-5f;
        }
        memcpy(aec_inst->post_res_frame + frame_shift, res_float,
               frame_shift * sizeof(float));
        memcpy(aec_inst->post_echo_frame + frame_shift, echo_float,
               frame_shift * sizeof(float));

        AecResidualEchoCancellation(
            aec_inst->post_res_frame, aec_inst->post_echo_frame,
            aec_inst->post_fft_conf, aec_inst->coh_res, aec_inst->coh_echo,
            aec_inst->coh_rd, aec_inst->adpt_coh_echo,
            aec_inst->prev_enchanced_sqrd, aec_inst->enhanced,
            &aec_inst->first_frame);
        memmove(aec_inst->post_res_frame,
                aec_inst->post_res_frame + frame_shift,
                sizeof(float) * frame_shift);
        memmove(aec_inst->post_echo_frame,
                aec_inst->post_echo_frame + frame_shift,
                sizeof(float) * frame_shift);
        for (q = 0; q < frame_shift; q += 2) {
            tmpno1 = aec_inst->vlp_buf[q] + aec_inst->enhanced[q];
            tmpno2 = aec_inst->vlp_buf[q + 1] + aec_inst->enhanced[q + 1];
            out[q] = (short)(tmpno1 * gain);
            out[q + 1] = (short)(tmpno2 * gain);
            aec_inst->vlp_buf[q] = aec_inst->enhanced[POST_FFT_LEN / 2 + q];
            aec_inst->vlp_buf[q + 1] =
                aec_inst->enhanced[POST_FFT_LEN / 2 + q + 1];
        }

        fwrite(out, sizeof(short), frame_shift, fid_out);
    }

    aec_inst->AudioProcessing_AEC_Release();

    fclose(fid_res);
    fclose(fid_out);
    fclose(fid_echo);

    if (CloseDelayRecordFile() < 0) {
        return -1;
    }

    return 0;
}

#ifdef __cplusplus
}
#endif
