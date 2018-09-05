#include "utility/fft/cdft_8g.h"
#include "utility/fft/bitrv2.h"
#include "utility/fft/bitrv2conj.h"
#include "utility/fft/ft_ops.h"

void cdft(int n, int isgn, float *a, int *ip, float *w) {
    if (n > 4) {
        if (isgn >= 0) {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
        } else {
            bitrv2conj(n, ip + 2, a);
            cftbsub(n, a, w);
        }
    } else if (n == 4) {
        cftfsub(n, a, w);
    }
}
