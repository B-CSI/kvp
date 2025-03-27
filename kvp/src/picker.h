#ifndef KVP_PK_H
#define KVP_PK_H

#include <stdint.h>
#include <stdbool.h>

typedef struct {
     int64_t idx_ref;
     int64_t idx_max;
     int64_t idx_ons;
     double t_ons;
     double posix_ons;
     double dt_jump;
     double kv_jump;
     void *KVPick; /* Placeholder for some wizzard-type capture refinement that will likely never happen */
} t_Singlepick;

typedef struct {
    t_Singlepick **picks;
    int *idx_bands;
    int *nbands;
    int nbands_max;
    int N_max;
} t_KVPicks;

int captures (t_KVPicks *grouped, t_Singlepick **refined, int *npicks, double *gaps, bool **mask);

#endif