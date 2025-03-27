#ifndef KVP_WT_H
#define KVP_WT_H

#define	PI 3.14159265358979328
#define NSIGMAS 8

typedef struct {
	double **wr;
	double *scale;           /* Vector of scales.             */
	unsigned int *Ls;        /* The length of each wavelet.   */
	int *center;             /* The wavelet origin.           */
	unsigned int *Down_smp;  /* Downsampling of each wavelet. */
	double Cpsi;             /* Normalization constant.       */
	double op1;              /* Mexhat parameter.             */
    unsigned int Ns;         /* Length of vector of scales.   */
} t_WaveletFamily;

typedef struct {
	double **d;
	unsigned int *N;
	unsigned int S;
} t_RWTvar;

void setscales0 (double *ps, double s0, unsigned int S, unsigned int V, double a0);
unsigned int setwaveletlength0 (unsigned int *Ls, double *ps, unsigned int S, unsigned int MaxL);
void setsampling0 (unsigned int *DownSmpl, unsigned int S, int continuous, unsigned int V, unsigned int J, double b0, double a0, double s0);

void MexicanHatFun_real (double *w, int *center, double dump, unsigned int Ls, double scale);
void FillWaveletFamily (t_WaveletFamily *pWF);

void cdotx_dd (double *y, double *x0, unsigned int N, double *w0, unsigned int L, int c, int inc);
int real_1D_wavelet_dec (t_RWTvar *y, double *x, unsigned int N, t_WaveletFamily *pWF);

#endif