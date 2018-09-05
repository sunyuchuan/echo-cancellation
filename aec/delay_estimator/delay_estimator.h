#ifndef DelayEstimator_MODULES_AUDIO_PROCESSING_UTILITY_DELAY_ESTIMATOR_H_
#define DelayEstimator_MODULES_AUDIO_PROCESSING_UTILITY_DELAY_ESTIMATOR_H_

typedef struct {
    // Pointer to bit counts.
    int* far_bit_counts;
    // Binary history variables.
    unsigned int* binary_far_history;
    int history_size;
} BinaryDelayEstimatorFarend;

typedef struct {
    // Pointer to bit counts.
    unsigned int* mean_bit_counts;
    // Array only used locally in ProcessBinarySpectrum() but whose size is
    // determined at run-time.
    unsigned int* bit_counts;

    // Binary history variables.
    unsigned int* binary_near_history;
    int near_history_size;

    // Delay estimation variables.
    unsigned int minimum_probability;
    int last_delay_probability;

    // Delay memory.
    int last_delay;

    // Far-end binary spectrum history buffer etc.
    BinaryDelayEstimatorFarend* farend;
} BinaryDelayEstimator;

int OpenDelayRecordFile();
int CloseDelayRecordFile();

// Releases the memory allocated by
// DelayEstimator_CreateBinaryDelayEstimatorFarend(...).
// Input:
//    - self              : Pointer to the binary delay estimation far-end
//                          instance which is the return value of
//                          DelayEstimator_CreateBinaryDelayEstimatorFarend().
//
void DelayEstimator_FreeBinaryDelayEstimatorFarend(
    BinaryDelayEstimatorFarend* self);

// Allocates the memory needed by the far-end part of the binary delay
// estimation. The memory needs to be initialized separately through
// DelayEstimator_InitBinaryDelayEstimatorFarend(...).
//
// Inputs:
//      - history_size    : Size of the far-end binary spectrum history.
//
// Return value:
//      - BinaryDelayEstimatorFarend*
//                        : Created |handle|. If the memory can't be allocated
//                          or if any of the input parameters are invalid NULL
//                          is returned.
//
BinaryDelayEstimatorFarend* DelayEstimator_CreateBinaryDelayEstimatorFarend(
    int history_size);

// Initializes the delay estimation far-end instance created with
// DelayEstimator_CreateBinaryDelayEstimatorFarend(...).
//
// Input:
//    - self              : Pointer to the delay estimation far-end instance.
//
// Output:
//    - self              : Initialized far-end instance.
//
void DelayEstimator_InitBinaryDelayEstimatorFarend(
    BinaryDelayEstimatorFarend* self);

// Adds the binary far-end spectrum to the internal far-end history buffer. This
// spectrum is used as reference when calculating the delay using
// DelayEstimator_ProcessBinarySpectrum().
//
// Inputs:
//    - self                  : Pointer to the delay estimation far-end
//                              instance.
//    - binary_far_spectrum   : Far-end binary spectrum.
//
// Output:
//    - self                  : Updated far-end instance.
//
void DelayEstimator_AddBinaryFarSpectrum(BinaryDelayEstimatorFarend* self,
                                         unsigned int binary_far_spectrum);

// Releases the memory allocated by
// DelayEstimator_CreateBinaryDelayEstimator(...).
//
// Note that BinaryDelayEstimator utilizes BinaryDelayEstimatorFarend, but does
// not take ownership of it, hence the BinaryDelayEstimator has to be torn down
// before the far-end.
//
// Input:
//    - self              : Pointer to the binary delay estimation instance
//                          which is the return value of
//                          DelayEstimator_CreateBinaryDelayEstimator().
//
void DelayEstimator_FreeBinaryDelayEstimator(BinaryDelayEstimator* self);

// Allocates the memory needed by the binary delay estimation. The memory needs
// to be initialized separately through
// DelayEstimator_InitBinaryDelayEstimator(...).
//
// Inputs:
//      - farend        : Pointer to the far-end part of the Binary Delay
//                        Estimator. This memory has to be created separately
//                        prior to this call using
//                        DelayEstimator_CreateBinaryDelayEstimatorFarend().
//
//                        Note that BinaryDelayEstimator does not take
//                        ownership of |farend|.
//
//      - lookahead     : Amount of non-causal lookahead to use. This can
//                        detect cases in which a near-end signal occurs before
//                        the corresponding far-end signal. It will delay the
//                        estimate for the current block by an equal amount,
//                        and the returned values will be offset by it.
//
//                        A value of zero is the typical no-lookahead case.
//                        This also represents the minimum delay which can be
//                        estimated.
//
//                        Note that the effective range of delay estimates is
//                        [-|lookahead|,... ,|history_size|-|lookahead|)
//                        where |history_size| was set upon creating the far-end
//                        history buffer size.
//
// Return value:
//      - BinaryDelayEstimator*
//                        : Created |handle|. If the memory can't be allocated
//                          or if any of the input parameters are invalid NULL
//                          is returned.
//
BinaryDelayEstimator* DelayEstimator_CreateBinaryDelayEstimator(
    BinaryDelayEstimatorFarend* farend, int lookahead);

// Initializes the delay estimation instance created with
// DelayEstimator_CreateBinaryDelayEstimator(...).
//
// Input:
//    - self              : Pointer to the delay estimation instance.
//
// Output:
//    - self              : Initialized instance.
//
void DelayEstimator_InitBinaryDelayEstimator(BinaryDelayEstimator* self);

// Estimates and returns the delay between the binary far-end and binary near-
// end spectra. It is assumed the binary far-end spectrum has been added using
// DelayEstimator_AddBinaryFarSpectrum() prior to this call. The value will be
// offset by the lookahead (i.e. the lookahead should be subtracted from the
// returned value).
//
// Inputs:
//    - self                  : Pointer to the delay estimation instance.
//    - binary_near_spectrum  : Near-end binary spectrum of the current block.
//
// Output:
//    - self                  : Updated instance.
//
// Return value:
//    - delay                 :  >= 0 - Calculated delay value.
//                              -2    - Insufficient data for estimation.
//
int DelayEstimator_ProcessBinarySpectrum(BinaryDelayEstimator* self,
                                         unsigned int binary_near_spectrum);

// Returns the last calculated delay updated by the function
// DelayEstimator_ProcessBinarySpectrum(...).
//
// Input:
//    - self                  : Pointer to the delay estimation instance.
//
// Return value:
//    - delay                 :  >= 0 - Last calculated delay value
//                              -2    - Insufficient data for estimation.
//
int DelayEstimator_binary_last_delay(BinaryDelayEstimator* self);

// Returns the estimation quality of the last calculated delay updated by the
// function DelayEstimator_ProcessBinarySpectrum(...). The estimation quality is
// a value in the interval [0, 1] in Q14. The higher the value, the better
// quality.
//
// Input:
//    - self                  : Pointer to the delay estimation instance.
//
// Return value:
//    - delay_quality         :  >= 0 - Estimation quality (in Q14) of last
//                                      calculated delay value.
//                              -2    - Insufficient data for estimation.
//
int DelayEstimator_binary_last_delay_quality(BinaryDelayEstimator* self);

#endif
