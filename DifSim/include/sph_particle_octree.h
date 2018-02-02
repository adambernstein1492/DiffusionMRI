#ifndef SPH_PARTICLE_OCTREE
#define SPH_PARTICLE_OCTREE
#include "sph_point.h"
#include "sph_particle_list.h"

#define SPH_PARTICLE_OCTREE_SPLIT_THRESH 16

struct _sph_particle_octree;
struct _sph_particle_octree_node;

typedef struct _sph_particle_octree_parent
{
  int total_particles;
  struct _sph_particle_octree_node *children;
} sph_particle_octree_parent;
typedef struct _sph_particle_octree_node
{
  union {
    sph_particle_list leaf;
    sph_particle_octree_parent parent;
  };
} sph_particle_octree_node;

typedef struct _sph_particle_octree
{
  float region_ranges[6]; // x,y,z min, x,y,z max
  sph_particle_octree_node root;
  sph_particle_data *dat;
} sph_particle_octree;

bool check_octree(sph_particle_octree *oct, int maxval);
void print_octree(sph_particle_octree *oct);
void init_octree(sph_particle_octree_node *root);
bool add_particle(sph_particle_octree *oct, sph_particle *p);
bool add_particle(sph_particle_octree_node *root, sph_particle *p,
                  double *p_offset, double *oct_width, sph_particle_data *dat);
void clear_octree(sph_particle_octree *oct);
#endif


