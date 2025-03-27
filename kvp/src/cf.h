#ifndef KVP_CF_H
#define KVP_CF_H

#define THRESHOLD_NAIVE 15.0

typedef struct {
	double **x;
	double **y;
	unsigned int *dims;
	unsigned int *L;
	unsigned int N;
} t_Movarray;

void trace_movkurtosis (double *y, double *const x, const unsigned int N, const unsigned int L);
void trace_movmxjump (double *const jmp, double *const x, const unsigned int N, const unsigned int L);

void movkurtosis (t_Movarray *cf_set);
void movjumps (t_Movarray *jmp_set);

#endif