#include <math.h>
#include <stdlib.h>

#include "filter.h"


/* Set the values of the scales sampled. */
void setscales0 (double *ps, double s0, unsigned int S, unsigned int V, double a0) {
	double da1, da2;
	unsigned int s;
	
	da1 = s0;
	da2 = pow(a0, 1/(double)V);
	for (s=0; s<S; s++) {
		ps[s] = da1;
		da1 *= da2;
	}
}

/* Set wavelet functions length. */
unsigned int setwaveletlength0 (unsigned int *Ls, double *ps, unsigned int S, unsigned int MaxL) {
	unsigned int s, Nsmpl, ua1;
	
	Nsmpl = 0;
	for (s=0; s<S; s++) {
		ua1 = 2*(unsigned int)ceil(NSIGMAS*ps[s]) + 1;
		Ls[s] = (ua1>MaxL) ? MaxL : ua1;
		Nsmpl += ua1;
	}
	return Nsmpl;
}

/* Set Downsampling field. */
void setsampling0 (unsigned int *DownSmpl, unsigned int S, int continuous, unsigned int V, unsigned int J, double b0, double a0, double s0) {
	double da1;
	unsigned int s, j, v, ua1;
	
	if (continuous) 
		for (s=0; s<S; s++) DownSmpl[s] = 1;
	else {
		da1 = s0*b0;
		for (j=0; j<J; j++) {
			ua1 = (da1 < 1) ? 1 : (unsigned int)da1;
			for (v=0; v<V; v++) *DownSmpl++ = ua1;
			da1 *= a0;
		}
	}
}

void MexicanHatFun_real (double *w, int *center, double dump, unsigned int Ls, double scale) {
	double k, daux;
	int u, u0;

	k = 2/sqrt(3*sqrt(PI)*scale);
	*center = Ls/2;
	scale = 1/scale;
	u0 = Ls/2;
	for (u = -u0; u < (signed)Ls - u0; u++) {
		daux  = scale*u;
		daux *= daux;
		*w++  = k *(daux - 1)*exp(-daux/2);
	}
}

/* Fill wavelet family frames. */
void FillWaveletFamily (t_WaveletFamily *pWF) {
	
	unsigned int s;
	
	for (s=1; s<pWF->Ns; s++) pWF->wr[s] = pWF->wr[s-1] + pWF->Ls[s-1];
	for (s=0; s<pWF->Ns; s++) {
		MexicanHatFun_real(pWF->wr[s], pWF->center + s, pWF->op1, pWF->Ls[s], pWF->scale[s]);
	}

	pWF->Cpsi = ( (4. / 3.) * sqrt(PI));
}

void cdotx_dd (double *y, double *x0, unsigned int N, double *w0, unsigned int L, int c, int inc) {
	double *x, *w, aux;
	unsigned int l, l0, n, n0;

	if (N<L) {
		x = x0; x0 = w0; w0 = x;
		l = N; N = L; L = l;
	}
	for (n=0; n<N; n+=inc) {
		n0 = (N + n - c) % N;
		l0 = N - n0;
		if (L < l0) l0 = L;

		x = x0 + n0;
		w = w0;
		aux = 0;
		for (l=0; l<l0>>1; l++) {
			aux += w[0]*x[0] + w[1]*x[1];
			w+=2; x+=2;
		}
		if (l0&1) aux += *w++*x[0];
		x = x0;
		for (l=0; l<(L-l0)>>1; l++) {
			aux += w[0]*x[0] + w[1]*x[1];
			w+=2; x+=2;
		}
		if ((L-l0)&1) aux += w[0]*x[0];
		*y++ = aux;
	}
}

int real_1D_wavelet_dec (t_RWTvar *y, double *x, unsigned int N, t_WaveletFamily *pWF) {
	unsigned int s, L;
	int er = 0;

	if (x == NULL) return -1;
	
	y->S = pWF->Ns;
	#pragma omp parallel for default(shared) private(L) schedule(dynamic,1)
	for (s=0; s<pWF->Ns; s++) {
		L = pWF->Ls[s];
		if (N < L) {
			er = -3; 
			continue;
		}
		cdotx_dd(y->d[s], x, N, pWF->wr[s], L, pWF->center[s], pWF->Down_smp[s]);
		y->N[s] = (N + pWF->Down_smp[s] - 1)/pWF->Down_smp[s];
	}
	return er;
}
