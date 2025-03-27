#include "cf.h"
#include <omp.h>


void movkurtosis (t_Movarray *cf_set) {
    // Moving window kurtosis for an arbitrary  number of traces
    
    double *x, *y;
    unsigned int size, L, N, i;
    
    N = cf_set->N;
    
    #pragma omp parallel for private(x, y, size, L)
    for (i=0; i<N; i++) {
        y = cf_set->y[i];
        x = cf_set->x[i];
        size = cf_set->dims[i];
        L = cf_set->L[i];
        trace_movkurtosis(y, x, size, L);
    }
    
}


void movjumps (t_Movarray *jmp_set) {
    // Maximum jump within a sliding window for an arbitrary number of traces
    
    double *x, *jmp;
    unsigned int size, L, N, i;
    
    N = jmp_set->N;
    
    // Parallelization happens inside "trace_movmxjump"
    for (i=0; i<N; i++) {
        jmp = jmp_set->y[i];
        x = jmp_set->x[i];
        size = jmp_set->dims[i];
        L = jmp_set->L[i];
        trace_movmxjump(jmp, x, size, L);
    }
    
}
