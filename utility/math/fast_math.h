#ifndef _FAST_MATH_H_
#define _FAST_MATH_H_

#define fmax(a,b) (((a)>(b))?(a):(b))
#define fmin(a,b) (((a)<(b))?(a):(b))

static __inline float _reciprocal_sqrt(float x)
{
	const float threeHalves = 1.5f;
	float number,xhalf,xcopy,tmp1,tmp2,tmp3,tmp4,tmp5;
	int xbitcopy;

	number = x;
	xhalf = 0.5f*x;
	xbitcopy = *(int*)&number;
	xbitcopy = 0x5f3759df - ( xbitcopy>>1 );
	xcopy   = * ( float * ) &xbitcopy;
	tmp1 = xcopy*xcopy;
	tmp2 = xcopy*threeHalves;
	tmp3 = xcopy*xhalf;
	tmp4 = tmp2 - tmp1*tmp3;
	tmp5 = tmp4*tmp4;
	tmp2 = tmp4*threeHalves;
	tmp3 = tmp4*xhalf;
	
	return (tmp2-tmp5*tmp3);
}

static __inline float _reciprocal_sqrt_hp(float x)
{
	const float threeHalves = 1.5f;
	float number,xhalf,xcopy,tmp1,tmp2,tmp3,tmp4,tmp5;
	int xbitcopy;

	number = x;
	xhalf = 0.5f*x;
	xbitcopy = *(int*)&number;
	xbitcopy = 0x5f3759df - ( xbitcopy>>1 );
	xcopy   = * ( float * ) &xbitcopy;
	tmp1 = xcopy*xcopy;
	tmp2 = xcopy*threeHalves;
	tmp3 = xcopy*xhalf;
	tmp4 = tmp2 - tmp1*tmp3;
	tmp5 = tmp4*tmp4;
	tmp2 = tmp4*threeHalves;
	tmp3 = tmp4*xhalf;

	tmp1 = tmp2-tmp5*tmp3;
	tmp2 = tmp1*tmp1;
	tmp3 = tmp1*threeHalves;
	tmp4 = tmp1*xhalf;

	return (tmp3-tmp2*tmp4);
}

static __inline float _reciprocal(float x)
{
	const float threeHalves = 1.5f;
	float number,xhalf,xcopy,tmp1,tmp2,tmp3,tmp4,tmp5;
	int xbitcopy;

	number = x*x;
	xhalf = 0.5f*number;
	xbitcopy = *(int*)&number;
	xbitcopy = 0x5f3759df - ( xbitcopy>>1 );
	xcopy   = * ( float * ) &xbitcopy;
	tmp1 = xcopy*xcopy;
	tmp2 = xcopy*threeHalves;
	tmp3 = xcopy*xhalf;
	tmp4 = tmp2 - tmp1*tmp3;
	tmp5 = tmp4*tmp4;
	tmp2 = tmp4*threeHalves;
	tmp3 = tmp4*xhalf;

	return (tmp2-tmp5*tmp3);
}

static __inline float _reciprocal_hp(float x)
{
	const float threeHalves = 1.5f;
	float number,xhalf,xcopy,tmp1,tmp2,tmp3,tmp4,tmp5;
	int xbitcopy;

	number = x*x;
	xhalf = 0.5f*number;
	xbitcopy = *(int*)&number;
	xbitcopy = 0x5f3759df - ( xbitcopy>>1 );
	xcopy   = * ( float * ) &xbitcopy;
	tmp1 = xcopy*xcopy;
	tmp2 = xcopy*threeHalves;
	tmp3 = xcopy*xhalf;
	tmp4 = tmp2 - tmp1*tmp3;
	tmp5 = tmp4*tmp4;
	tmp2 = tmp4*threeHalves;
	tmp3 = tmp4*xhalf;

	tmp1 = tmp2-tmp5*tmp3;
	tmp2 = tmp1*tmp1;
	tmp3 = tmp1*threeHalves;
	tmp4 = tmp1*xhalf;

	return (tmp3-tmp2*tmp4);
}

static __inline float _ln(const float val, register float* const pTable, register const unsigned precision)
{
	/* get access to float bits */
	register const int* const pVal = (const int*)(&val);

	/* extract exponent and mantissa (quantized) */
	register const int exp = ((*pVal >> 23) & 255) - 127;
	register const int man = (*pVal & 0x7FFFFF) >> (23 - precision);

	/* exponent plus lookup refinement */
	return ((float)(exp) + pTable[man]) * 0.693147180559945f;
}

static __inline float _exp(const float val, unsigned int* pTable,const unsigned int  precision)
{
	// build float bits
	const int i = (int) (val * 12102203.16156145696768f + 1065353216.0f);

	// replace mantissa with lookup
	const int it = (i & 0xFF800000) | pTable[(i & 0x7FFFFF) >> (23 - precision)];

	// convert bits to float
	union { int i; float f; } pun;
	return pun.i = it,  pun.f;
}

static __inline float _tanh(const float val,register unsigned int* pTable,const unsigned int  precision)
{
	register const float exp_val = _exp(2.0f*val, pTable, precision);
	//register const float recp_exp = _reciprocal(exp_val+1.0f);
	register const float exp_val2_minus1 = exp_val - 1.0f;
	register const float exp_val2_plus1 = exp_val + 1.0f;
	return exp_val2_minus1/exp_val2_plus1;
	//return recp_exp*exp_val2_minus1;
}


#endif
