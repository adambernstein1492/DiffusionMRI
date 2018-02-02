#ifndef SPH_PARTICLE
#define SPH_PARTICLE
#include "sph_point.h"
typedef struct _sph_particle
{
  sph_point loc, lastloc, orig_loc;
  float mass, prev_density, prev_press;
  char flags;
} sph_particle_data;

// now just an index
typedef int sph_particle;
#endif
