// Performs delay estimation on block by block basis.
// The return value is  0 - OK and -1 - Error, unless otherwise stated.

#ifndef _AUDIO_PROCESSING_UTILITY_DELAY_ESTIMATOR_WRAPPER_H_
#define _AUDIO_PROCESSING_UTILITY_DELAY_ESTIMATOR_WRAPPER_H_

// Releases the memory allocated by
// DelayEstimator_CreateDelayEstimatorFarend(...) Input:
//      - handle        : Pointer to the delay estimation far-end instance.
//
void DelayEstimator_FreeDelayEstimatorFarend(void* handle);

// Allocates the memory needed by the far-end part of the delay estimation. The
// memory needs to be initialized separately through
// DelayEstimator_InitDelayEstimatorFarend(...).
//
// Inputs:
//      - spectrum_size : Size of the spectrum used both in far-end and
//                        near-end. Used to allocate memory for spectrum
//                        specific buffers.
//      - history_size  : The far-end history buffer size. Note that the maximum
//                        delay which can be estimated is controlled together
//                        with |lookahead| through
//                        DelayEstimator_CreateDelayEstimator().
//
// Return value:
//      - void*         : Created |handle|. If the memory can't be allocated or
//                        if any of the input parameters are invalid NULL is
//                        returned.
//
void* DelayEstimator_CreateDelayEstimatorFarend(int spectrum_size,
                                                int history_size);

// Initializes the far-end part of the delay estimation instance returned by
// DelayEstimator_CreateDelayEstimatorFarend(...)
// Input:
//      - handle        : Pointer to the delay estimation far-end instance.
//
// Output:
//      - handle        : Initialized instance.
//
int DelayEstimator_InitDelayEstimatorFarend(void* handle);

// Adds the far-end spectrum to the far-end history buffer. This spectrum is
// used as reference when calculating the delay using
// DelayEstimator_ProcessSpectrum().
//
// Inputs:
//    - handle          : Pointer to the delay estimation far-end instance.
//    - far_spectrum    : Far-end spectrum.
//    - spectrum_size   : The size of the data arrays (same for both far- and
//                        near-end).
//    - far_q           : The Q-domain of the far-end data.
//
// Output:
//    - handle          : Updated far-end instance.
//
int DelayEstimator_AddFarSpectrumFix(void* handle, float* far_spectrum,
                                     int spectrum_size);

// See DelayEstimator_AddFarSpectrumFix() for description.
int DelayEstimator_AddFarSpectrumFloat(void* handle, float* far_spectrum,
                                       int spectrum_size, int band_first,
                                       int band_last);

// Releases the memory allocated by DelayEstimator_CreateDelayEstimator(...)
// Input:
//      - handle        : Pointer to the delay estimation instance.
//
void DelayEstimator_FreeDelayEstimator(void* handle);

// Allocates the memory needed by the delay estimation. The memory needs to be
// initialized separately through DelayEstimator_InitDelayEstimator(...).
//
// Inputs:
//      - farend_handle : Pointer to the far-end part of the delay estimation
//                        instance created prior to this call using
//                        DelayEstimator_CreateDelayEstimatorFarend().
//
//                        Note that DelayEstimator_CreateDelayEstimator does not
//                        take ownership of |farend_handle|, which has to be
//                        torn down properly after this instance.
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
//      - void*         : Created |handle|. If the memory can't be allocated or
//                        if any of the input parameters are invalid NULL is
//                        returned.
//
void* DelayEstimator_CreateDelayEstimator(void* farend_handle, int lookahead);

// Initializes the delay estimation instance returned by
// DelayEstimator_CreateDelayEstimator(...)
// Input:
//      - handle        : Pointer to the delay estimation instance.
//
// Output:
//      - handle        : Initialized instance.
//
int DelayEstimator_InitDelayEstimator(void* handle);

// Estimates and returns the delay between the far-end and near-end blocks. The
// value will be offset by the lookahead (i.e. the lookahead should be
// subtracted from the returned value).
// Inputs:
//      - handle        : Pointer to the delay estimation instance.
//      - near_spectrum : Pointer to the near-end spectrum data of the current
//                        block.
//      - spectrum_size : The size of the data arrays (same for both far- and
//                        near-end).
//      - band_first    : From this band to start estimating delay
//      - band_last     : To this band to stop estimating delay
//
// Output:
//      - handle        : Updated instance.
//
// Return value:
//      - delay         :  >= 0 - Calculated delay value.
//                        -1    - Error.
//                        -2    - Insufficient data for estimation.
//
int DelayEstimator_DelayEstimatorProcessFloat(void* handle,
                                              float* near_spectrum,
                                              int spectrum_size, int band_first,
                                              int band_last);

// Returns the last calculated delay updated by the function
// DelayEstimator_DelayEstimatorProcess(...).
//
// Input:
//      - handle        : Pointer to the delay estimation instance.
//
// Return value:
//      - delay         : >= 0  - Last calculated delay value.
//                        -1    - Error.
//                        -2    - Insufficient data for estimation.
//
int DelayEstimator_last_delay(void* handle);

// Returns the estimation quality/probability of the last calculated delay
// updated by the function DelayEstimator_DelayEstimatorProcess(...). The
// estimation quality is a value in the interval [0, 1] in Q9. The higher the
// value, the better quality.
//
// Input:
//      - handle        : Pointer to the delay estimation instance.
//
// Return value:
//      - delay_quality : >= 0  - Estimation quality (in Q9) of last calculated
//                                delay value.
//                        -1    - Error.
//                        -2    - Insufficient data for estimation.
//
int DelayEstimator_last_delay_quality(void* handle);

#endif
