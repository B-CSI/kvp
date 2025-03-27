#include "cf.h"


inline void moments (double *M4, double *M3, double *M2, double *x, double mn, const unsigned int N) {
	double da1, da2, m2=0, m3=0, m4=0;
	unsigned int n;
	
	for (n=0; n<N; n++) {
		da1 = x[n] - mn;
		da2 = da1*da1;
		m2 += da2;
		m3 += da2*da1;
		m4 += da2*da2;
	}
	*M2 = m2;
	*M3 = m3;
	*M4 = m4;
}

void trace_movkurtosis (double *y, double *const x, const unsigned int N, const unsigned int L) {
	double da0, da1, M2=0, M3=0, M4=0, invL, invLp1;
	unsigned int n, L2 = L>>1, Lw;
	double m, m1, m2, cm, dm;
	
	/* Initialization and first sample */
	cm = 0;
	for (n=0; n<L2+1; n++) cm += x[n];
	m = cm/(double)n;
	moments (&M4, &M3, &M2, x, m, L2+1);
	*y++ = n*M4/(M2*M2);
	
	/* Starting samples: the window is partially outside the sequence domain. */
	for (   ; n<L; n++) {
		m   = cm/(double)n;
		dm  = (n*x[n] - cm)/(double)(n*(n+1));
		cm += x[n];
		m1  = cm/(double)(n+1);
		da0 = (x[n] - m1)*(x[n] - m);
		M4 += dm*(-4*M3 + 6*M2*dm + da0*dm*(n*n - n + 1));
		M3 += dm*(-3*M2 + da0*(n-1));
		M2 += da0;
		
		da0 = (n+1)*M4/(M2*M2);
		if (da0 > 15 || da0 < 1.5) moments (&M4, &M3, &M2, x, m, L);
		*y++ = (n+1)*M4/(M2*M2);
	}
	
	/* Central samples: the window is inside the sequence domain. */
	invL = 1/(double)L;
	invLp1 = 1/(double)(L+1);
	for (   ; n<N; n++) {
		m1  = cm*invL;
		m2  = (cm  + x[n])*invLp1;
		dm  = (L*x[n] - cm)*invLp1*invL; 
		cm += x[n] - x[n-L];
		m   = cm*invL;
		da0 = (x[n] - m1)*(x[n] - m2);
		da1 = (x[n-L] - m2)*(x[n-L] - m);
		M4 += dm*(-4*M3 + 6*M2*dm + da0*dm*(L*L - L + 1));
		M3 += dm*(-3*M2 + da0*(L-1));
		M2 += da0 - da1;
		dm = (L*x[n-L] - cm)*invLp1*invL;
		M3 -= dm*(-3*M2 + da1*(L-1));
		M4 -= dm*(-4*M3 + 6*M2*dm + da1*dm*(L*L - L + 1));
		
		da0 = L*M4/(M2*M2);
		if (da0 > THRESHOLD_NAIVE || da0 < 1.5) moments (&M4, &M3, &M2, x+n-L+1, m, L);
		*y++ = L*M4/(M2*M2);
	}
	/* Ending samples. */
	for (   ; n<N+L2; n++) {
		Lw = L+N-n-1;
		m1   = cm/(Lw+1);
		cm  -= x[n-L];
		m    = cm/Lw;
		da0 = (x[n-L] - m1)*(x[n-L] - m);
		M2 -= da0;
		M3 -= -3*M2*(m1-m) + (m1-m)*da0*(Lw-1);
		M4 -= -4*M3*(m1-m) + 6*M2*(m1-m)*(m1-m) + da0*(m1-m)*(m1-m)*(Lw*Lw - Lw + 1);
		*y++ = Lw*M4/(M2*M2);
	}
	
}
