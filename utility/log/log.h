#ifndef _AEC_LOG_H_
#define _AEC_LOG_H_

#ifdef __ANDROID__
#include <android/log.h>

#define TAG "debug-AEC"
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, TAG, __VA_ARGS__)
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO, TAG, __VA_ARGS__)
#define LOGW(...) __android_log_print(ANDROID_LOG_WARN, TAG, __VA_ARGS__)
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR, TAG, __VA_ARGS__)
#define LOGF(...) __android_log_print(ANDROID_LOG_FATAL, TAG, __VA_ARGS__)
#elif __linux__
#define LOGD printf
#define LOGI printf
#define LOGW printf
#define LOGE printf
#define LOGF printf
#endif

#endif
