#include <omp.h>


inline void minmax (unsigned int *imn, unsigned int *imx, double * const x, const unsigned int B, const unsigned int E) {
	unsigned int l, mn, mx;
	double da1;
	
	mn = B; mx = B;
	for (l=B+1; l<=E; l++) {
		da1 = x[l];
		if (x[mn] > da1) mn = l;
		else if (x[mx] < da1) mx = l;
	}
	*imn = mn; *imx = mx;
}

inline unsigned int max (double * const x, const unsigned int B, const unsigned int E) {
	unsigned int l, mx;
	
	mx = B;
	for (l=B+1; l<=E; l++)
		if (x[mx] < x[l]) mx = l;
		
	return mx;
}

inline unsigned int min (double * const x, const unsigned int B, const unsigned int E) {
	unsigned int l, mn;
	
	mn = B;
	for (l=B+1; l<=E; l++)
		if (x[mn] > x[l]) mn = l;
		
	return mn;
}

void trace_movmxjump (double *const jmp, double *const x, const unsigned int N, const unsigned int L) {
	unsigned int n, L2=(L+1)/2, imx, imn, B, E;
	double da1;
	
	minmax (&imn, &imx, x, 0, L2-1);
	jmp[0] = (imn > imx) ? x[imn]-x[imx]:x[imx]-x[imn];
	
	for (n=1; n<L-L2; n++) {
		E = n+L2-1;
		
		da1 = x[E];
		if (x[imn] > da1) imn = E;
		else if (x[imx] < da1) imx = E;
		
		jmp[n] = (imn > imx) ? x[imn]-x[imx]:x[imx]-x[imn];
	}
#if 0   /* Sequential version */
	for (n=L-L2; n<N-(L-L2); n++) {
		B = n-(L-L2);
		E = B+L-1;
		
		if (imx < B) imx = max (x, B, E);
		else if (x[imx] < x[E]) imx = E;
		
		if (imn < B) imn = min (x, B, E);
		else if (x[imn] > x[E]) imn = E;
		
		jmp[n] = (imn > imx) ? x[imn]-x[imx]:x[imx]-x[imn];
	}
#else   /* Parallel version: the key issues are the dependences across loop iterations (mx & mn) */
	#pragma omp parallel default(shared) private(B, E, n)
	{
		unsigned int n1, n2, Step, mn, mx;
		unsigned int Nth = omp_get_num_threads();
		unsigned int Idth = omp_get_thread_num();
		unsigned int Elements = N-2*(L-L2)+1;
		
		Step = (Elements+Nth-1)/Nth;
		n1 = Step*Idth + L-L2;
		n2 = n1 + Step;
		n2 = (n2>Elements + L-L2) ? Elements  + L-L2 : n2;
		
		B = n1-(L-L2);
		E = B+L-1;
		minmax (&mn, &mx, x, B, E);
		
		for (n=n1; n<n2; n++) {
			B = n-(L-L2);
			E = B+L-1;
			
			if (mx < B) mx = max (x, B, E);
			else if (x[mx] < x[E]) mx = E;
			
			if (mn < B) mn = min (x, B, E);
			else if (x[mn] > x[E]) mn = E;
			
			jmp[n] = (mn > mx) ? x[mn]-x[mx]:x[mx]-x[mn];
		}
		if (Idth == Nth-1) { imn = mn; imx = mx; }
	}
#endif
	for (n=N-(L-L2); n<N; n++) {
		B = n-(L-L2);
		E = N-1;
		
		if (imx < B) imx = max (x, B, E);
		else if (imn < B) imn = min (x, B, E);
		
		jmp[n] = (imn > imx) ? x[imn]-x[imx]:x[imx]-x[imn];
	}
	
}