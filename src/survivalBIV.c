#include <math.h>
#include <R.h>
#include <string.h>

/* FIXME:
  Code extracted from https://cran.r-project.org/src/contrib/Archive/condSURV/ - 2022-04-29
*/

/*
Author:
	Artur Araújo (adapted from 'rcmp')

Description:
	Compares two double values.

Parameters:
	x[in]		first double
	y[in]		second double

Return value:
	Returns -1 if x is lower than y.
	Returns 1 if x is greater than y.
	Returns 0 otherwise.
*/

static int cmp_doubles(
	const double x,
	const double y)
{
	if (x < y) return -1;
	if (x > y) return 1;
	return 0;
} // cmp_doubles


/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Sorts 'x' by increasing order with 'indx', 'y' and 'indy' alongside.

Parameters:
	x[inout]		pointer to vector's first element
	indx[inout]		pointer to vector's first element
	y[inout]		pointer to vector's first element
	indy[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_biv(
	double *const x,
	int *const indx,
	double *const y,
	int *const indy,
	const int n)
{
	double v, u;
	int i, j, h, iv, iu;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			v = x[i];
			iv = indx[i];
			u = y[i];
			iu = indy[i];
			j = i;
			while (j >= h && cmp_doubles(x[j - h], v) > 0) {
				x[j] = x[j-h];
				indx[j] = indx[j-h];
				y[j] = y[j-h];
				indy[j] = indy[j-h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
			y[j] = u;
			indy[j] = iu;
		}
	}
} // sort_biv

/*
Author:
	Artur Araújo (adapted from 'icmp')

Description:
	Compares two int values.

Parameters:
	x[in]		first int
	y[in]		second int

Return value:
	Returns -1 if x is lower than y.
	Returns 1 if x is greater than y.
	Returns 0 otherwise.
*/

static int cmp_ints(
	const int x,
	const int y)
{
	if (x < y) return -1;
	if (x > y) return 1;
	return 0;
} // cmp_ints

/*
Author:
	Artur Araújo (adapted from 'rsort_with_index')

Description:
	Restores the elements of a sorted double vector back
		to their previous order.

Parameters:
	indx[inout]		pointer to index vector's first element
	x[inout]		pointer to vector's first element
	n[in]			vector's length

Return value:
	This function doesn't return a value.
*/

static void sort_back_double(
	int *const indx,
	double *const x,
	const int n)
{
	double v;
	int i, j, h, iv;

	for (h = 1; h <= n / 9; h = 3 * h + 1);
	for (; h > 0; h /= 3) {
		for (i = h; i < n; i++) {
			iv = indx[i];
			v = x[i];
			j = i;
			while (j >= h && cmp_ints(indx[j - h], iv) > 0) {
				indx[j] = indx[j-h];
				x[j] = x[j-h];
				j -= h;
			}
			indx[j] = iv;
			x[j] = v;
		}
	}
} // sort_back_double

/*
Author:
	Artur Araújo

Description:
	Computes weights based on a kernel.

Parameters:
	time1[in]		pointer to covariate first element
	len[in]			pointer to length of covariate
	t1[in]			pointer to covariate value to compute the weights at
	lambda[in]		pointer to lambda
	kernel[in]		pointer to pointer that points to a char array
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.
*/

void WeightsKernel(
	const double *const time1,
	const int *const len,
	const double *const t1,
	const double *const lambda,
	const char *const *const kernel,
	double *const weights)
{
	char *kernel1 = "gaussian";
	char *kernel2 = "epanechnikov";
	char *kernel3 = "tricube"; // also known as triweight
	char *kernel4 = "boxcar"; // also known as uniform
	char *kernel5 = "triangular";
	char *kernel6 = "quartic"; // also known as biweight
	char *kernel7 = "cosine";
	register int i;
	for (i = 0; i < *len; i++) {
		weights[i] = (time1[i]-*t1) / *lambda;
	}
	if (strcmp(*kernel, kernel1) == 0) { // if kernel is gaussian
		for (i = 0; i < *len; i++) {
			weights[i] = exp(-pow(weights[i], 2)/2);
		}
	}
	else if (strcmp(*kernel, kernel2) == 0) { // if kernel is epanechnikov
		for (i = 0; i < *len; i++) {
			weights[i] = ( 1-pow(weights[i], 2) )*(fabs(weights[i]) <= 1);
		}
	}
	else if (strcmp(*kernel, kernel3) == 0) { // if kernel is tricube
		for (i = 0; i < *len; i++) {
			weights[i] = fabs(weights[i]);
			weights[i] = pow(1-pow(weights[i], 3), 3)*(weights[i] <= 1);
		}
	}
	else if (strcmp(*kernel, kernel4) == 0) { // if kernel is boxcar
		for (i = 0; i < *len; i++) {
			weights[i] = (fabs(weights[i]) <= 1);
		}
	}
	else if (strcmp(*kernel, kernel5) == 0) { // if kernel is triangular
		for (i = 0; i < *len; i++) {
			weights[i] = fabs(weights[i]);
			weights[i] = (1-weights[i])*(weights[i] <= 1);
		}
	}
	else if (strcmp(*kernel, kernel6) == 0) { // if kernel is quartic
		for (i = 0; i < *len; i++) {
			weights[i] = pow(1-pow(weights[i], 2), 2)*(fabs(weights[i]) <= 1);
		}
	}
	else if (strcmp(*kernel, kernel7) == 0) { // if kernel is cosine
		for (i = 0; i < *len; i++) {
			weights[i] = cos(M_PI*weights[i]/2)*(fabs(weights[i]) <= 1);
		}
	}
	return;
} // WeightsKernel

/*
Author:
	Artur Araújo

Description:
	Computes nearest neighbour weights.

Parameters:
	time1[in]		pointer to covariate first element
	len[in]			pointer to length of covariate
	t1[in]			pointer to covariate value to compute the weights at
	span[in]		pointer to span
	window[in]		pointer to pointer that points to a char array
	weights[out]	pointer to weights first element

Return value:
	This function doesn't return a value.
*/

void WeightsNNE(
	const double *const time1,
	const int *const len,
	const double *const t1,
	const double *const span,
	const char *const *const window,
	double *const weights)
{
	char *window1 = "asymmetric";
	char *window2 = "symmetric";
	register int i;
	int index[*len], tr;
	double lambda1;
	for (i = 0; i < *len; i++) {
		weights[i] = time1[i]-*t1; // initialize weights vector
		index[i] = i; // initialize index vector
	}
	rsort_with_index(weights, index, *len); // use internal R sorting to sort weights
	if (strcmp(*window, window1) == 0) { // if window is asymmetric
		tr = *len * *span;
		for (i = *len-1; i > 0; i--) {
			if (weights[i] < 0) break; // compute index
		}
		if (i+tr+1 > *len-1) lambda1 = weights[*len-1];
		else lambda1 = weights[i+tr+1];
		for (i = 0; i < *len; i++) {
			weights[i] = (weights[i] >= 0 && weights[i] <= lambda1); // compute weights
		}
	}
	else if (strcmp(*window, window2) == 0) { // if window is symmetric
		double lambda0;
		tr = *len * *span / 2;
		for (i = *len-1; i > 0; i--) {
			if (weights[i] <= 0) break; // compute lower index
		}
		if (i-tr < 0) lambda0 = -fabs(weights[0]);
		else lambda0 = -fabs(weights[i-tr]);
		for (; i > 0; i--) {
			if (weights[i] < 0) break; // compute upper index
		}
		if (i+tr+1 > *len-1) lambda1 = weights[*len-1];
		else lambda1 = weights[i+tr+1];
		for (i = 0; i < *len; i++) {
			weights[i] = ( (weights[i] >= lambda0 && weights[i] <= 0) || (weights[i] >= 0 && weights[i] <= lambda1) ); // compute weights
		}
	}
	sort_back_double(index, weights, *len); // put weights back to their previous order
	return;
} // WeightsNNE

/*
Author:
	Artur Araújo

Description:
	Computes a single probability value at a specified time index.
	This function implements the weights version of the Kaplan-Meier
		estimator.

Parameters:
	time2[in]		pointer to time2 first element
	status[in]		pointer to status first element
	weights[in]		pointer to weights first element
	delta[in]		pointer to delta first element
	len[in]			pointer to length of time2, status, weights, and delta
	end[in]			pointer to index to compute the probability at
	surv[out]		pointer to survival probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, weights and delta must be sorted first.
	Vectors time2, status, weights and delta must have the same length.
*/

void WeightedKaplanMeierValue(
	const double *const time2,
	const int *const status,
	const double *const weights,
	const int *const delta,
	const int *const len,
	const int *const end,
	double *const surv)
{
	register int i;
	double n, d;
	*surv = 1;
	for (i = *len-1, n = 0; i >= *end; i--) { // loop in reverse order until end index is reached
		n += delta[i]*weights[i]; // initialize living weighting
	}
	while (i >= 0) { // loop through the sample in reverse order until zero index is reached
		n += delta[i]*weights[i];
		d = status[i]*weights[i]; // initialize dead weighting
		for (i--; i >= 0 && time2[i] == time2[i+1]; i--) { // loop in reverse order until time changes or zero index is reached
			n += delta[i]*weights[i]; // weight the living
			d += status[i]*weights[i]; // weight the dead
		}
		if (n > 0) *surv *= 1-d/n; // compute survival probability
	}
	return;
} // WeightedKaplanMeierValue

/*
Author:
	Artur Araújo

Description:
	Sorts time2, status, weights and delta and then calls 'WeightedKaplanMeierValue'.

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	weights[inout]	pointer to weights first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time2, status, weights, and delta
	t[in]			pointer to time to compute the probability at (defines last index)
	surv[out]		pointer to survival probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, weights and delta must have the same length.
*/

void WeightedKaplanMeierValueSort(
	double *const time2,
	int *const status,
	double *const weights,
	int *const delta,
	const int *const len,
	const double *const t,
	double *const surv)
{
	int end = *len/2;
	sort_biv(time2, status, weights, delta, *len); // sort data
	if (time2[end] > *t) end = 0;
	for (; end < *len; end++) {
		if (time2[end] > *t) break; // determine last index
	}
	WeightedKaplanMeierValue(time2, status, weights, delta, len, &end, surv); // compute survival probability
	return;
} // WeightedKaplanMeierValueSort

/*
Author:
	Artur Araújo

Description:
	Computes the conditional survival probability P(T2>t2|T1=t1).

Parameters:
	time2[inout]	pointer to time2 first element
	status[inout]	pointer to status first element
	time1[inout]	pointer to time1 first element
	delta[inout]	pointer to delta first element
	len[in]			pointer to length of time2, status, time1 and delta
	t2[in]			pointer to time2 value to compute the probability at
	t1[in]			pointer to time1 value to compute the probability at
	lambda[in]		pointer to lambda
	kernel[in]		pointer to pointer that points to a char array
	surv[out]		pointer to survival probability value

Return value:
	This function doesn't return a value.

Remarks:
	Vectors time2, status, time1 and delta must have the same length.
*/

void SurvBeranKernel(
	double *const time2,
	int *const status,
	double *const time1,
	int *const delta,
	const int *const len,
	const double *const t2,
	const double *const t1,
	const double *const lambda,
	const char *const *const kernel,
	double *const surv)
{
	double weights[*len];
	WeightsKernel(time1, len, t1, lambda, kernel, weights); // compute kernel weights
	WeightedKaplanMeierValueSort(time2, status, weights, delta, len, t2, surv); // compute conditional survival probability
	return;
} // SurvBeranKernel