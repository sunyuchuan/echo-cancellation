#ifndef _AEC_LOG_H_
#define _AEC_LOG_H_

#ifdef DEBUG_IN_PC

#define LOGD printf
#define LOGI printf
#define LOGW printf
#define LOGE printf
#define LOGF printf

#else
#include <android/log.h>

#define TAG "debug-AEC"
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, TAG, __VA_ARGS__)
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO, TAG, __VA_ARGS__)
#define LOGW(...) __android_log_print(ANDROID_LOG_WARN, TAG, __VA_ARGS__)
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR, TAG, __VA_ARGS__)
#define LOGF(...) __android_log_print(ANDROID_LOG_FATAL, TAG, __VA_ARGS__)
#endif

#endif
