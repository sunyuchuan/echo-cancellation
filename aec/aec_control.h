#ifndef _AEC_CONTROL_H_
#define _AEC_CONTROL_H_

#include "aec/post_process/post_process_config.h"
#include "subband/analy_synth/filterbank_control.h"
#include "utility/fft/fft_config.h"
#include "utility/ringbuffer/ring_buffer.h"

namespace xmly_audio_recorder_android {
class AecControl {
   public:
    int far_buf_write_pos;
    int far_buf_read_pos;
    int known_delay;
    int last_known_delay;

    short *far_frame_buf;
    short *near_frame_buf;

    RingBuffer *far_end_buf;
    short *last_far_end_frame;

    void *delay_estimator_farend;
    void *delay_estimator;
    int band_first;
    int band_last;
    unsigned short current_delay;
    float *far_history;
    float *far_sepctrum_history;
    int far_history_pos;

    FFT_Config *realFFT;
    FilterBankControl *fb_ctrl_far;
    FilterBankControl *fb_ctrl_near;
    FilterBankControl *fb_ctrl_echo;

    float *near_input_buf;
    float *far_input_buf;

    float *adp_filter_coeff;
    float *far_fb_buf;
    float *near_fb_buf;
    float *echo_buf;
    float *error_buf;
    float *post_res_frame;
#if AEC_POST_PROCESSING_ON
    FFT_Config *post_fft_conf;
#if AEC_POST_PROCESSING_COH
    float nonlinearity;
    float *coh_res;
    float *coh_echo;
    float *coh_rd;
    float *adpt_coh_echo;
    float *prev_enchanced_sqrd;
#endif
#if AEC_POST_PROCESSING_NN
    float *ln_lookup;
    int   *exp_lookup;
    float *joint_buf;
    float *nn_buf;
#endif
    float *enhanced;
    short first_frame;
    float *vlp_buf;

#endif
    float *post_echo_frame;
#if	AEC_SAVE_FAR
    float *far_frame_save;
#endif

    float *Rss;
    float *Rdd;
    float *Ree;

    AecControl();
    /**************************************************************************
     *Function  -   allocate space for aec instance
     *
     *Input		-	None
     *Output	-	None
     *
     *Return	-	-1 - fail
     *				0  - succeed
     **************************************************************************/
    int AudioProcessing_AEC_Create();

    /**************************************************************************
     *Function  -   initialize aec instance data members
     *
     *Input	    -	amp_perc: percentage of amplifier(0.0f ~ 1.0f)
     *          -	min_perc: min percentage of amplifier(0.0f ~ 1.0f)
     *          -	delay_band_select - 0 - use low band
     *		    -	1 - use high band
     *Output	-	None
     *
     *Return	-	-1 - fail
     *				0  - succeed
     **************************************************************************/
    int AudioProcessing_AEC_Init(float amp_perc, float min_perc,
                                 int delay_band_select);

    /**************************************************************************
     *Function  -   fill aec far-end buffer
     *
     *Input		-	ref - far-end data addr
     *				sample_size - far-end data length in byte
     *				playout_switch - audio playout switch
     *				mic_switch - recorder microphone switch
     *Output	-	None
     *
     *Return	-	-1 - fail
     *				0  - succeed
     **************************************************************************/
    int AudioProcessing_AEC_FillFarBuf(char *ref, short sample_size,
                                       bool playout_switch, bool mic_switch);

    /**************************************************************************
     *Function  -   adaptive filter process
     *
     *Input		-	pri - near-end data addr
     *			sample_size - near-end data length in byte
     *			out_buf - output data addr
     *			output_size - output data size addr in byte
     *			out_buf1 - estimated echo data addr
     *			output_size1 - estimated echo data size addr in
     *			playout_switch - if playout thread is working
     *			mic_switch - if microphone is capturing
     *			amp_pwr_changed - if the amplifier power is manually
     *changed amp_level - current amplifier power level min_level - from what
     *level to add nonlinearity byte playout_switch - audio playout switch
     *mic_switch - microphone capture switch Output	-	None
     *
     *Return	-	-1 - fail
     *				0  - succeed
     **************************************************************************/
    int AudioProcessing_AEC_Process(char *pri, short sample_size, char *out_buf,
                                    unsigned int *output_size, char *out_buf1,
                                    unsigned int *output_size1,
                                    bool playout_switch, bool mic_switch,
                                    bool amp_pwr_changed, float amp_level,
                                    float min_level);

    /**************************************************************************
     *Function  -   reset far frame buffer variables
     *
     *Input		-	None
     *Output	-	None
     *
     *Return	-	-1 - fail
     *				0  - succeed
     **************************************************************************/
    int AudioProcessing_AEC_ClearFarFrameBuf();

    /**************************************************************************
     *Function  -   reset near end filter bank state
     *
     *Input		-	None
     *Output	-	None
     *
     *Return	-	-1 - fail
     *				0  - succeed
     **************************************************************************/
    int AudioProcessing_AEC_ResetNearState();

    /**************************************************************************
     *Function  -   reset far end buffer and filter bank state
     *
     *Input		-	None
     *Output	-	None
     *
     *Return	-	-1 - fail
     *				0  - succeed
     **************************************************************************/
    int AudioProcessing_AEC_ResetFarState();

    /**************************************************************************
     *Function  -   release aec instance
     *
     *Input		-	None
     *Output	-	None
     *
     *Return	-	-1 - fail
     *				0  - succeed
     **************************************************************************/
    int AudioProcessing_AEC_Release();
};
}  // namespace xmly_audio_recorder_android

#endif
