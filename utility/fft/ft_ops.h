#ifndef _FOURIER_TRANSFORM_OPERATIONS_H_
#define _FOURIER_TRANSFORM_OPERATIONS_H_

void cftfsub(int n, float *a, float *w);
void cftbsub(int n, float *a, float *w);

void rftfsub(int n, float *a, int nc, float *c);
void rftbsub(int n, float *a, int nc, float *c);

#endif