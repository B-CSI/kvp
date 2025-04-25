#include "picker.h"
#include <math.h>


/* The compiler complains if these two functions are inlined, no idea why */
void wipe_auxiliaries (t_Singlepick **aux_p, int *aux_i, int n) {
    
    int i;
    
    for (i=0; i<n; i++) {
        aux_p[i] = 0; // NULL is not a constant on gcc 9.4 (Ubuntu 20.04), so 0 it is
        aux_i[i] = 0;
    }
    
}

void wipe_idx_buffers (int *idx, int nb) {
    
    int i;
    
    for (i=0; i<nb; i++) idx[i] = 0;
    
}

int captures (t_KVPicks *grouped, t_Singlepick **refined, int *npicks, double *gaps, bool **mask) {
    /* An unused Python version of this procedure (see callers.py) fully describes this algorithm */
	/* This algorithm removes a big bottleneck caused by the former, extremely naive one */
	
    int nb_max=grouped->nbands_max;
    
    t_Singlepick **picks=grouped->picks;
    int *idx_bands=grouped->idx_bands, *nbands=grouped->nbands;
    
    t_Singlepick *pick_ptrs_aux[nb_max];
    int idx_bands_aux[nb_max];
    
    int idx_p=0, idx_aux[nb_max];
    int idx_b=0, idx_r=0, idx_c=0, idx_s=0, idx_m=0;
    
    t_Singlepick *ref_band, *cas_band;
    double gap;
    bool *ref_mask, *cas_mask;
    
    int nb, search_from, cap_idx=0;
    double ref_ons, dt, dt_abs;
    
    int i=0;
    
    wipe_auxiliaries (pick_ptrs_aux, idx_bands_aux, nb_max);
    wipe_idx_buffers (idx_aux, nb_max);
    
    for (idx_b=0; idx_b<nb_max; idx_b++) {
        
        ref_band = refined[idx_b];
        ref_mask = mask[idx_b];
        gap = gaps[idx_b];
        
        for (idx_r=0; idx_r<npicks[idx_b]; idx_r++) {
            
            if (ref_mask[idx_r]) continue;
            
            nb = 1;
            
            idx_bands_aux[0] = idx_b;
            pick_ptrs_aux[0] = &ref_band[idx_r];
            ref_mask[idx_r] = 1;
            
            ref_ons = ref_band[idx_r].t_ons;
            
            for (idx_c=idx_b+1; idx_c<nb_max; idx_c++) {
                
                cas_band = refined[idx_c];
                cas_mask = mask[idx_c];
                search_from = idx_aux[idx_c];
                
                cap_idx = -1;
                
                for (idx_s=search_from; idx_s<npicks[idx_c]; idx_s++) {
                    
                    if (cas_mask[idx_s]) continue;
                    
                    dt = ref_ons - cas_band[idx_s].t_ons;
                    
                    if (dt > gap) continue;
                    
                    if (dt < -gap) {
                        idx_aux[idx_c] = idx_s;
                        break;
                    }
                    
                    /* If the iteration makes it here, there is a capture */
                    
                    dt_abs = fabs(dt);
                    cap_idx = idx_s;
                    
                    for (idx_m=idx_s+1 ;npicks[idx_c]; idx_m++) {
                        
                        if (cas_mask[idx_m]) break;
                        
                        dt = fabs(ref_ons-cas_band[idx_m].t_ons);
                        
                        if (dt > dt_abs) break;
                        
                        cap_idx = idx_m;
                        
                    }
                    
                    idx_bands_aux[nb] = idx_c;
                    pick_ptrs_aux[nb] = &cas_band[cap_idx];
                    cas_mask[cap_idx] = 1;
                    nb++;
                    
                    idx_aux[idx_c] = cap_idx;
                    
                    break;
                    
                }
                
            }
            
            nbands[idx_p] = nb;
            
            for (i=0; i<nb_max; i++) {
                picks[i] = pick_ptrs_aux[i];
                idx_bands[i] = idx_bands_aux[i];
            }
            
            wipe_auxiliaries (pick_ptrs_aux, idx_bands_aux, nb_max);
            
            idx_p++;
            picks += nb_max;
            idx_bands += nb_max;
            
        }
        
        wipe_idx_buffers (idx_aux, nb_max);
          
    }
    
    return idx_p;
      
}
