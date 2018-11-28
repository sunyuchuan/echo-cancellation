#include "aec/delay_estimator/delay_estimator_wrapper.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "aec/delay_estimator/delay_estimator.h"
#include "aec/delay_estimator/delay_estimator_internal.h"
//#include "compile_assert.h"

// Only bit |kBandFirst| through bit |kBandLast| are processed and
// |kBandFirst| - |kBandLast| must be < 32.
// enum { kBandFirst = 12 };
// enum { kBandLast = 43 };

// enum { kBandFirst = 7 };  // for 44k
// enum { kBandLast = 38 };

// enum { kBandFirst = 96 };  // for 44k
// enum { kBandLast = 127 };

float ar_factor1 = 0.0625f;

static __inline unsigned int SetBit(unsigned int in, int pos) {
    unsigned int mask = (1 << pos);
    unsigned int out = (in | mask);

    return out;
}

static __inline void DelayEstimator_MeanEstimatorFloat(float new_value,
                                                       float factor,
                                                       float* mean_value) {
    float diff = new_value - *mean_value;

    // mean_new = mean_value + ((new_value - mean_value) * factor);

    *mean_value += diff * factor;
}

// Computes the binary spectrum by comparing the input |spectrum| with a
// |threshold_spectrum|. Float and fixed point versions.
//
// Inputs:
//      - spectrum            : Spectrum of which the binary spectrum should be
//                              calculated.
//      - threshold_spectrum  : Threshold spectrum with which the input
//                              spectrum is compared.
// Return:
//      - out                 : Binary spectrum.
//
static unsigned int BinarySpectrumFix(float* spectrum,
                                      SpectrumType* threshold_spectrum,
                                      int* threshold_initialized,
                                      int kBandFirst, int kBandLast) {
    int i = kBandFirst;
    unsigned int out = 0;

    if (!(*threshold_initialized)) {
        // Set the |threshold_spectrum| to half the input |spectrum| as starting
        // value. This speeds up the convergence.
        for (i = kBandFirst; i <= kBandLast; i++) {
            if (spectrum[i] > 0) {
                threshold_spectrum[i].float_ = spectrum[i] * 0.5f;
                *threshold_initialized = 1;
            }
        }
    }
    for (i = kBandFirst; i <= kBandLast; i++) {
        // Update the |threshold_spectrum|.
        DelayEstimator_MeanEstimatorFloat(spectrum[i], ar_factor1,
                                          &(threshold_spectrum[i].float_));
        // Convert |spectrum| at current frequency bin to a binary value.
        if (spectrum[i] > threshold_spectrum[i].float_) {
            out = SetBit(out, i - kBandFirst);
        }
    }

    return out;
}

void DelayEstimator_FreeDelayEstimatorFarend(void* handle) {
    DelayEstimatorFarend* self = (DelayEstimatorFarend*)handle;

    if (handle == NULL) {
        return;
    }

    free(self->mean_far_spectrum);
    self->mean_far_spectrum = NULL;

    DelayEstimator_FreeBinaryDelayEstimatorFarend(self->binary_farend);
    self->binary_farend = NULL;

    free(self);
}

void* DelayEstimator_CreateDelayEstimatorFarend(int spectrum_size,
                                                int history_size) {
    DelayEstimatorFarend* self = NULL;

    // Check if the sub band used in the delay estimation is small enough to fit
    // the binary spectra in a uint32_t.
    // COMPILE_ASSERT(kBandLast - kBandFirst < 32);

    // if (spectrum_size >= kBandLast) {
    //    self = malloc(sizeof(DelayEstimator));
    //}
    self = malloc(sizeof(DelayEstimator));

    if (self != NULL) {
        int memory_fail = 0;

        // Allocate memory for the binary far-end spectrum handling.
        self->binary_farend =
            DelayEstimator_CreateBinaryDelayEstimatorFarend(history_size);
        memory_fail |= (self->binary_farend == NULL);

        // Allocate memory for spectrum buffers.
        self->mean_far_spectrum = malloc(spectrum_size * sizeof(SpectrumType));
        memory_fail |= (self->mean_far_spectrum == NULL);

        self->spectrum_size = spectrum_size;

        if (memory_fail) {
            DelayEstimator_FreeDelayEstimatorFarend(self);
            self = NULL;
        }
    }

    return self;
}

int DelayEstimator_InitDelayEstimatorFarend(void* handle) {
    DelayEstimatorFarend* self = (DelayEstimatorFarend*)handle;

    if (self == NULL) {
        return -1;
    }

    // Initialize far-end part of binary delay estimator.
    DelayEstimator_InitBinaryDelayEstimatorFarend(self->binary_farend);

    // Set averaged far and near end spectra to zero.
    memset(self->mean_far_spectrum, 0,
           sizeof(SpectrumType) * self->spectrum_size);
    // Reset initialization indicators.
    self->far_spectrum_initialized = 0;

    return 0;
}

int DelayEstimator_AddFarSpectrumFloat(void* handle, float* far_spectrum,
                                       int spectrum_size, int band_first,
                                       int band_last) {
    DelayEstimatorFarend* self = (DelayEstimatorFarend*)handle;
    unsigned int binary_spectrum = 0;

    if (self == NULL) {
        return -1;
    }
    if (far_spectrum == NULL) {
        // Empty far end spectrum.
        return -1;
    }
    if (spectrum_size != self->spectrum_size) {
        // Data sizes don't match.
        return -1;
    }

    // Get binary spectrum.
    binary_spectrum = BinarySpectrumFix(far_spectrum, self->mean_far_spectrum,
                                        &(self->far_spectrum_initialized),
                                        band_first, band_last);
    DelayEstimator_AddBinaryFarSpectrum(self->binary_farend, binary_spectrum);

    return 0;
}

void DelayEstimator_FreeDelayEstimator(void* handle) {
    DelayEstimator* self = (DelayEstimator*)handle;

    if (handle == NULL) {
        return;
    }

    free(self->mean_near_spectrum);
    self->mean_near_spectrum = NULL;

    DelayEstimator_FreeBinaryDelayEstimator(self->binary_handle);
    self->binary_handle = NULL;

    free(self);
}

void* DelayEstimator_CreateDelayEstimator(void* farend_handle, int lookahead) {
    DelayEstimator* self = NULL;
    DelayEstimatorFarend* farend = (DelayEstimatorFarend*)farend_handle;

    if (farend_handle != NULL) {
        self = malloc(sizeof(DelayEstimator));
    }

    if (self != NULL) {
        int memory_fail = 0;

        // Allocate memory for the farend spectrum handling.
        self->binary_handle = DelayEstimator_CreateBinaryDelayEstimator(
            farend->binary_farend, lookahead);
        memory_fail |= (self->binary_handle == NULL);

        // Allocate memory for spectrum buffers.
        self->mean_near_spectrum =
            malloc(farend->spectrum_size * sizeof(SpectrumType));
        memory_fail |= (self->mean_near_spectrum == NULL);

        self->spectrum_size = farend->spectrum_size;

        if (memory_fail) {
            DelayEstimator_FreeDelayEstimator(self);
            self = NULL;
        }
    }

    return self;
}

int DelayEstimator_InitDelayEstimator(void* handle) {
    DelayEstimator* self = (DelayEstimator*)handle;

    if (self == NULL) {
        return -1;
    }

    // Initialize binary delay estimator.
    DelayEstimator_InitBinaryDelayEstimator(self->binary_handle);

    // Set averaged far and near end spectra to zero.
    memset(self->mean_near_spectrum, 0,
           sizeof(SpectrumType) * self->spectrum_size);
    // Reset initialization indicators.
    self->near_spectrum_initialized = 0;

    return 0;
}

int DelayEstimator_DelayEstimatorProcessFloat(void* handle,
                                              float* near_spectrum,
                                              int spectrum_size, int band_first,
                                              int band_last) {
    DelayEstimator* self = (DelayEstimator*)handle;
    unsigned int binary_spectrum = 0;

    if (self == NULL) {
        return -1;
    }
    if (near_spectrum == NULL) {
        // Empty near end spectrum.
        return -1;
    }
    if (spectrum_size != self->spectrum_size) {
        // Data sizes don't match.
        return -1;
    }

    // Get binary spectra.
    binary_spectrum = BinarySpectrumFix(near_spectrum, self->mean_near_spectrum,
                                        &(self->near_spectrum_initialized),
                                        band_first, band_last);

    return DelayEstimator_ProcessBinarySpectrum(self->binary_handle,
                                                binary_spectrum);
}

int DelayEstimator_last_delay(void* handle) {
    DelayEstimator* self = (DelayEstimator*)handle;

    if (self == NULL) {
        return -1;
    }

    return DelayEstimator_binary_last_delay(self->binary_handle);
}

int DelayEstimator_last_delay_quality(void* handle) {
    DelayEstimator* self = (DelayEstimator*)handle;

    if (self == NULL) {
        return -1;
    }

    return DelayEstimator_binary_last_delay_quality(self->binary_handle);
}
