#include "utility/fft/rdft_8g_init.h"
#include "utility/fft/rdft_8g_make_table.h"

void rdft_init(int* ip, float* w, int n, int* nw, int* nc) {
    int nw1, nc1;
    nw1 = ip[0];

    if (n > (nw1 << 2)) {
        *nw = n >> 2;
        makewt(*nw, ip, w);
    }
    nc1 = ip[1];
    if (n > (nc1 << 2)) {
        *nc = n >> 2;
        makect(*nc, ip, w + (*nw));
    }
}
