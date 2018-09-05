#ifndef _RingBuffer_H_
#define _RingBuffer_H_

#include <stddef.h>  // size_t

enum Wrap { SAME_WRAP, DIFF_WRAP };

typedef struct RingBuffer {
    size_t read_pos;
    size_t write_pos;
    size_t element_count;
    size_t element_size;
    enum Wrap rw_wrap;
    char* data;
} RingBuffer;

// Returns NULL on failure.
RingBuffer* RingBuffer_CreateBuffer(size_t element_count, size_t element_size);
int RingBuffer_InitBuffer(RingBuffer* handle);
void RingBuffer_FreeBuffer(void* handle);

// Reads data from the buffer. The |data_ptr| will point to the address where
// it is located. If all |element_count| data are feasible to read without
// buffer wrap around |data_ptr| will point to the location in the buffer.
// Otherwise, the data will be copied to |data| (memory allocation done by the
// user) and |data_ptr| points to the address of |data|. |data_ptr| is only
// guaranteed to be valid until the next call to RingBuffer_WriteBuffer().
//
// To force a copying to |data|, pass a NULL |data_ptr|.
//
// Returns number of elements read.
size_t RingBuffer_ReadBuffer(RingBuffer* handle, void** data_ptr, void* data,
                             size_t element_count);

// Writes |data| to buffer and returns the number of elements written.
size_t RingBuffer_WriteBuffer(RingBuffer* handle, const void* data,
                              size_t element_count);

// Moves the buffer read position and returns the number of elements moved.
// Positive |element_count| moves the read position towards the write position,
// that is, flushing the buffer. Negative |element_count| moves the read
// position away from the the write position, that is, stuffing the buffer.
// Returns number of elements moved.
int RingBuffer_MoveReadPtr(RingBuffer* handle, int element_count);

// Returns number of available elements to read.
size_t RingBuffer_available_read(const RingBuffer* handle);

// Returns number of available elements for write.
size_t RingBuffer_available_write(const RingBuffer* handle);

#endif
