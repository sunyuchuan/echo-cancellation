#ifndef _AUDIO_PROCESSING_UTILITY_DELAY_ESTIMATOR_INTERNAL_H_
#define _AUDIO_PROCESSING_UTILITY_DELAY_ESTIMATOR_INTERNAL_H_

#include "aec/delay_estimator/delay_estimator.h"

typedef union {
    float float_;
    int int32_;
} SpectrumType;

typedef struct {
    // Pointers to mean values of spectrum.
    SpectrumType* mean_far_spectrum;
    // |mean_far_spectrum| initialization indicator.
    int far_spectrum_initialized;

    int spectrum_size;

    // Far-end part of binary spectrum based delay estimation.
    BinaryDelayEstimatorFarend* binary_farend;
} DelayEstimatorFarend;

typedef struct {
    // Pointers to mean values of spectrum.
    SpectrumType* mean_near_spectrum;
    // |mean_near_spectrum| initialization indicator.
    int near_spectrum_initialized;

    int spectrum_size;

    // Binary spectrum based delay estimator
    BinaryDelayEstimator* binary_handle;
} DelayEstimator;

#endif
