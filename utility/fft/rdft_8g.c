#include "utility/fft/rdft_8g.h"
#include "utility/fft/bitrv2.h"
#include "utility/fft/ft_ops.h"

void rdft(int n, int *nc, int *nw, int isgn, float *a, int *ip, float *w) {
    float xi;

    if (isgn >= 0) {
        if (n > 4) {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
            rftfsub(n, a, *nc, w + *nw);
        } else if (n == 4) {
            cftfsub(n, a, w);
        }
        xi = a[0] - a[1];
        a[0] += a[1];
        a[1] = xi;
    } else {
        a[1] = 0.5f * (a[0] - a[1]);
        a[0] -= a[1];
        if (n > 4) {
            rftbsub(n, a, *nc, w + *nw);
            bitrv2(n, ip + 2, a);
            cftbsub(n, a, w);
        } else if (n == 4) {
            cftfsub(n, a, w);
        }
    }
}
