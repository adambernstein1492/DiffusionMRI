
#ifndef SPH_PARTICLE_LIST
#define SPH_PARTICLE_LIST
#include "sph_particle.h"

typedef struct _sph_particle_list
{
  int nparticles, nalloc;
  sph_particle *particles;
} sph_particle_list;

void init_particle_list(sph_particle_list *pl);
bool rem_particle(sph_particle_list *plist, sph_particle *p);
bool add_particle(sph_particle_list *plist, sph_particle *p);
#endif
