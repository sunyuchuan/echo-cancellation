#include "utility/fft/rdft_8g_make_table.h"
#include <math.h>
#include "utility/fft/bitrv2.h"

void makewt(int nw, int *ip, float *w) {
    int j, nwh;
    float delta, x, y;

    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2) {
        nwh = nw >> 1;
        delta = atan(1.0f) / nwh;
        w[0] = 1;
        w[1] = 0;
        w[nwh] = cos(delta * nwh);
        w[nwh + 1] = w[nwh];
        if (nwh > 2) {
            for (j = 2; j < nwh; j += 2) {
                x = cos(delta * j);
                y = sin(delta * j);
                w[j] = x;
                w[j + 1] = y;
                w[nw - j] = y;
                w[nw - j + 1] = x;
            }
            for (j = nwh - 2; j >= 2; j -= 2) {
                x = w[2 * j];
                y = w[2 * j + 1];
                w[nwh + j] = x;
                w[nwh + j + 1] = y;
            }
            bitrv2(nw, ip + 2, w);
        }
    }
}

void makect(int nc, int *ip, float *c) {
    int j, nch;
    float delta;

    ip[1] = nc;
    if (nc > 1) {
        nch = nc >> 1;
        delta = atan(1.0f) / nch;
        c[0] = cos(delta * nch);
        c[nch] = 0.5f * c[0];
        for (j = 1; j < nch; j++) {
            c[j] = 0.5f * cos(delta * j);
            c[nc - j] = 0.5f * sin(delta * j);
        }
    }
}
