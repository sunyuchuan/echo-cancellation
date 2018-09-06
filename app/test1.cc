#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aec/aec_control.h"
#include "aec/aec_defines.h"
#include "aec/delay_estimator/delay_estimator.h"
#include "subband/analy_synth/filterbank.h"

#ifdef __cplusplus
}
#endif

using namespace xmly_audio_recorder_android;
int main() {
    AecControl *aec_inst = new AecControl;
    FILE *fid_ref, *fid_pri, *fid_out, *fid_echo;
    char *ref, *pri, out[4096], out1[4096];
    short block_size = 4096;
    unsigned int size, size1, count = 0;
    bool mic, playout;

    fid_ref = fopen("/home/layne/audio/record/test_aec/play.pcm", "rb");
    fid_pri = fopen("/home/layne/audio/record/test_aec/sdl.pcm", "rb");
    fid_out = fopen("out.pcm", "wb");
    fid_echo = fopen("echo.pcm", "wb");

    ref = (char *)malloc(block_size);
    pri = (char *)malloc(block_size);
    memset(ref, 0, block_size);
    memset(pri, 0, block_size);

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

    while (fread(ref, sizeof(char), block_size, fid_ref) == block_size &&
           fread(pri, sizeof(char), block_size, fid_pri) == block_size) {
        count++;
        aec_inst->AudioProcessing_AEC_FillFarBuf(ref, block_size, playout, mic);
        aec_inst->AudioProcessing_AEC_Process(pri, block_size, out, &size, out1,
                                              &size1, playout, mic);
        fwrite(out, sizeof(char), size, fid_out);
        fwrite(out1, sizeof(char), size1, fid_echo);
    }

    aec_inst->AudioProcessing_AEC_Release();

    fclose(fid_ref);
    fclose(fid_pri);
    fclose(fid_out);
    fclose(fid_echo);
    if (ref) {
        free(ref);
        ref = NULL;
    }
    if (pri) {
        free(pri);
        pri = NULL;
    }

    if (CloseDelayRecordFile() < 0) {
        return -1;
    }

    return 0;
}
