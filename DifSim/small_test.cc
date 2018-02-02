#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <getopt.h>
#include "ConfigFile.h"
#include "rng.h"

#define MAX_TESS_PTS 1
#define R_UINT_MAX 2.3283064365386963e-10
#define MY_PI 3.14159265358979323846

using std::ifstream;
using std::ofstream;
using std::ios;
using std::string;
using std::cout;
using std::endl;

struct input_data {
  int Delta;
  int delta;
  double gradient_strength;
  int ramp;
  int pulse;
  double snr;
  double voxel_size;
  double volume_size;
  double Dextra;
  int nsteps;
  int step_size;
  int density;
  int ntess;
  int fiber_facets;
  int number_of_particles;
  int simulation_type;
  string signal_file;
  string mcell_file;
};

double *r_step;
int radial_subdivisions = 1024;
struct rng_state *rng;

double r_func(double s)
{   
  double f,s_sqr,val;
    
  f=2.25675833419102511712;  /* 4.0/sqrt(pi) */
  s_sqr=s*s;
  val=f*s_sqr*exp(-s_sqr);
    
  return(val);
}   

int init_r_step(void)
{   
  double inc,target,accum,r,r_max,delta_r,delta_r2;
  int j;
    
  if ((r_step=(double *)malloc(radial_subdivisions*sizeof(double)))==NULL) { 
    cout << "MCell: cannot store radial step length table" << endl;
    return(1);
  } 
      
  inc=1.0/radial_subdivisions;
  accum=0;
  r_max=3.5;
  delta_r=r_max/(1000*radial_subdivisions);
  delta_r2=0.5*delta_r;
  r=0;
  target=0.5*inc;
  j=0;
  while (j<radial_subdivisions) {
   accum=accum+(delta_r*r_func(r+delta_r2));
   r=r+delta_r;
   if (accum>=target) {
     r_step[j]=r;
     target=target+inc;
     j++;
   }
  }
  printf("Min r step = %20.17g   Max r step = %20.17g\n",
            r_step[0],r_step[radial_subdivisions-1]);
  return(0);
}   

void parse_input(input_data *d, char *file_name) {
  ConfigFile config(file_name);
  d->Delta = config.read("Delta", 100000);
  d->delta = config.read("delta", 10000);
  d->gradient_strength = config.read("gradient strength", 4.0);
  d->ramp = config.read("ramp", 500);
  d->pulse = config.read("pulse", 0);
  d->snr = config.read("snr", 100);
  d->voxel_size = config.read("voxel size", 50.0);
  d->volume_size = config.read("volume size", 100.0);
  d->Dextra = config.read("Dextra", 0.075);
  d->nsteps = config.read("nsteps", 1);
  d->step_size = config.read("step size", 100);
  d->density = config.read("density", 16);
  d->ntess = config.read("ntess", 3);
  d->fiber_facets = config.read("fiber facets", 10);
  d->number_of_particles = config.read("nparts", 4096);
  d->simulation_type = config.read("simulation type", 1);
  if ( !config.readInto<std::string>(d->signal_file, "signal file",
				     "signal_file.txt") ) {
    cout << "Using default value for signal file" << endl;
  }
  config.readInto<std::string>(d->mcell_file, "mcell file",
			       "mcell.mdl");
}

void simulation(input_data *d) {
  int rank = 0;
  int n_machines = 1;

  double *particle_location;
  double *spin_signal; // summed over all directions, two dimensional
  double *vertex_signal;
  double *particle_phase_shifts; // for each particle and each gradient
  double time_unit;

  double voxel_size, volume_size, Dintra, Dextra, Tperma,  G;
  int NSTEPS, density, ntess, sph_slices, pulse, step_size, Delta, delta, ramp;
  double snr;
  int nparticles;
  int sim_type;

  Delta = d->Delta;
  delta = d->delta;
  G = d->gradient_strength;
  ramp = d->ramp;
  pulse = d->pulse;
  snr = d->snr;
  voxel_size = d->voxel_size;
  volume_size = d->volume_size;
  Dextra = d->Dextra;
  NSTEPS = d->nsteps;
  step_size = d->step_size;
  density = d->density;
  ntess = d->ntess;
  nparticles = d->number_of_particles;
  nparticles *= 1000;
  sim_type = d->simulation_type;

  int nvert;
  double *gx;
  double *gy;
  double *gz;
  int *tris;

  double r, h, r_sin_phi, theta;

  tris  = (int *) calloc((unsigned)(2*MAX_TESS_PTS),3*sizeof(int));

  particle_location = (double *)calloc(nparticles, sizeof(double));
  particle_phase_shifts = (double *)calloc(nparticles, sizeof(double));
  spin_signal = (double*)calloc(2*nvert, sizeof(double));
  vertex_signal = (double*)calloc(nvert, sizeof(double));
  
  double proj, factor, ang;
  if (ntess == 0) {
    nvert = 1;
  }

  init_r_step();

  /* time step and mc_factor for delta */
  double timestep = 1.0e-6;
  time_unit = delta*timestep;
  double mc_factor = sqrt(4.0*1.0e8*Dextra*time_unit);
  cout << "mc_factor: "  << mc_factor << endl;
  factor = G * 26752.0e-10;
  rng = (struct rng_state *)(malloc(sizeof(struct rng_state)));
  rng_init(rng, 1);

  for (int i = 0; i < nparticles; i++) {
    r = mc_factor*r_step[(int)(R_UINT_MAX*radial_subdivisions*rng_uint(rng))];
    h = (2.0*R_UINT_MAX*rng_uint(rng))-1.0;
    r_sin_phi=r*sqrt(1.0-(h*h));
    theta = 2*MY_PI*R_UINT_MAX*rng_uint(rng);
    particle_location[i]=r_sin_phi*cos(theta);
    proj = particle_location[i] * factor;
    particle_phase_shifts[i] += proj;
    ang = particle_phase_shifts[i];
    /*    if (i < 10) {
      cout << "ang " << i << ": " << ang << endl;
      cout << "r " << i << ": " << r << endl;
      }*/
    /*    spin_signal[0] += cos(ang);
	  spin_signal[1] += sin(ang); */
  }

  /* time step and mc_factor for Delta */
  time_unit = Delta*timestep;
  mc_factor = sqrt(4.0*1.0e8*Dextra*time_unit);
  cout << "mc_factor: "  << mc_factor << endl;

  double max_distance = 0.0;

  for (int i = 0; i < nparticles; i++) {
    r = mc_factor*r_step[(int)(R_UINT_MAX*radial_subdivisions*rng_uint(rng))];
    h = (2.0*R_UINT_MAX*rng_uint(rng))-1.0;
    r_sin_phi=r*sqrt(1.0-(h*h));
    theta = 2*MY_PI*R_UINT_MAX*rng_uint(rng);
    particle_location[i] += r_sin_phi*cos(theta);
    if (fabs(particle_location[i]) > max_distance) {
      max_distance = fabs(particle_location[i]);
    }
    proj = particle_location[i] * factor;
    particle_phase_shifts[i] += proj;
    ang = - particle_phase_shifts[i];
    /*
    if (i < 10) {
      cout << "ang " << i << ": " << ang << endl;
      cout << "r " << i << ": " << r << endl;
    }
    */
    spin_signal[0] += cos(ang);
    spin_signal[1] += sin(ang);
  }
  double renorm = 1.0/nparticles;
  vertex_signal[0] = hypot(spin_signal[0] * renorm,
			spin_signal[1] * renorm);
 
  cout << "spin signal: " << spin_signal[0] << ", " << spin_signal[1] << endl;
  cout << "vertex signal: " << vertex_signal[0] << endl;
  cout << "max diffusion: " << max_distance << endl;
}

static char *progname;
static const char *short_options = "h";
static const char *short_optind = "h";
static struct option long_options[] = {
  {"help",               no_argument, 0, 'h'},
  {(const char*)0, 0, 0, 0}
};

static void printUsage(std::ostream& out) {
  out << "Usage: " << progname << " Configuration_File " << endl;
}

static void parseOptions(int argc, char **argv)
{
  int option_index = 0;
  int next_option;
  while ((next_option = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1) {
  again:
    switch (next_option) {
    case 0:
      next_option = short_optind[option_index];
      goto again;
    case 'h': printUsage(std::cerr); exit(0);
      break;
    default:
      std::cerr << "Unknown option: " << (char)next_option << std::endl;
      printUsage(std::cerr);
      exit(1);
    }
  }
}

int main(int argc, char *argv[]) {
  input_data d;
  progname = argv[0];
  parseOptions(argc, argv);
  if (optind==argc) { printUsage(std::cerr); exit(0); }
  parse_input(&d, argv[1]);
  simulation(&d);
  return 0;
}
