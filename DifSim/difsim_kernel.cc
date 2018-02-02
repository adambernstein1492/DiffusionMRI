#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <numeric>
#include <time.h>

#include "mcell_wrap.h"
#include "mcell_structs.h"
#include "ConfigFile.h"
#include "gtb_histogram.h"
#include "VectorPotential.h"

/* Dot product between two vectors: x.y */
#define DOT_VEC3(x,y) ((x)[0]*(y)[0] + (x)[1]*(y)[1] + (x)[2]*(y)[2])
#define SUM_VEC3(a,x,y) a[0] = x[0] + y[0], a[1] = x[1] + y[1], a[2] = x[2] + y[2]
#define SCALE_VEC3(a,s,x) a[0] = s*x[0], a[1] = s*x[1], a[2] = s*x[2]

//extern struct ligand_info **ligand_table;
extern struct volume *world;

/* Pulse sequence types */
#define PULSE_GRAD_ECHO 0
#define PULSE_SPIN_ECHO 1
#define PULSE_DPFG 2
#define PULSE_DWSSFP 3
#define PULSE_MULTI_SPIN_ECHO 4
#define PULSE_STEAM_DTI 5

/* Sampling options */
static const char *SAMPLE_TYPES[] = \
  {"Polytope tessellation on sphere","Linear spherical","Random spherical"} ;

#define SAMP_TESS 0
#define SAMP_SPHL 1
#define SAMP_SRND 2

/* Tessellation polytopes */
static const char *POLY_TYPES[] = {"Octahedral","Tetrahedral","Icosahedral"} ;

#define POLY_OCT 0
#define POLY_TET 1
#define POLY_ICO 2

static int POLY_TYPE = POLY_ICO; /* default = icosahedral */

static char *progname;

#define MAX_TESS_PTS 20000

int sphtes( float *gx, float *gy, float *gz, int *tris,
            int polytope, int ntess );

using std::ostringstream;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::complex;
using std::stringstream;
using std::max;

typedef enum {NONE, HEX2D, CCP3D, HCP3D, ANYREGION} structure_t;

/* simulation variables */
struct mcellSim {
  int debugMCell;
  float timestep, orig_timestep;
  int num_particles;
  int gradientDirections;
  int gradient_on_time;
  int gradient_off_time;
  int gradient_ramp_time;
  int pulse_sequence; // 0 for spin echo, 1 for gradient echo
  int mixing_time;
  int endTime;
  int periodic;
  int threeRegion;
  int permeable;
  int parallelPlates;
  int cylinderArray;
  int structureArray;
  structure_t structurePacking;
  float bath_concentration;
  float sheath_concentration;
  float core_concentration;
  float D_bath;
  float D_sheath;
  float D_core;
  float current_gradient_strength;
  float gradient_strength; /* units G/cm */
  float voxel_size;
  float xDim;
  float yDim;
  float zDim;
  float xScale;
  float yScale;
  float zScale;
  float volume_size;
  float volumeHeight;
  float coreRadius;
  float sheathRadius;
  float sphereRadius;
  float cylinderRadius;
  float cylinderSpacing;
  float cylinderHeight;
  float dpfgAngle;
  dvec  dpfgAngleList;
  int   tensorFileDPFG;
  float EulerAngle1;
  float EulerAngle2;
  float EulerAngle3;
  float ncEulerAngle1;
  float ncEulerAngle2;
  float ncEulerAngle3;
  float P_scale;
  float ncCurrent;
  float Tr[3][3];
  float ncTr[3][3];
  int ncNx;
  int ncNy;
  int ncNz;
  CVectorPotential *vp;
  dvec  ncCurrentValues;
  dvec  ncCurrentTimes;
  int repeat_number;
  int repeat_time;
  int flip_angle;
  int use_relaxation;
  int echo_time;
  int echo_spacing;
  int echo_number;
  bc_t BCx0;
  bc_t BCx1;
  bc_t BCy0;
  bc_t BCy1;
  bc_t BCz0;
  bc_t BCz1;
  dvec gx1;
  dvec gy1;
  dvec gz1;
  dvec gx2;
  dvec gy2;
  dvec gz2;
};

/* pulse sequence timings
   
   these times represent the leading edge of each phase of the pulse sequnce */
struct sequenceTimes {
  int rampUp1, rampUp2;
  int rampDown1, rampDown2;
  int gOn1, gOn2;
  int gOff;
  int end;
};

/* input variables */
struct input_data {
  int debugMCell;
  int Delta;
  int delta;
  int ramp;
  int pulse;
  int mixing_time;
  int nsteps;
  double step_size;
  int ntess;
  int gradientDirections;
  int number_of_particles;
  int randomSeed;
  int threeRegion;
  int permeable;
  int cylinderArray;
  int structureArray;
  structure_t structurePacking;
  string signal_filename;
  string complex_signal_filename;
  string mcell_file;
  string log_file;
  string tensor_file;
  string structureFile;
  string sheathFile;
  string ncStrengthFile;
  string ncOutputFile;
  string ncObjectName;
  string ncJinFile;
  string ncJoutFile;
  string dpfgAngleFile;
  float voxel_size;
  float xDim;
  float yDim;
  float zDim;
  float xScale;
  float yScale;
  float zScale;
  float volume_size;
  float volumeHeight;
  float coreRadius;
  float sheathRadius;
  float sphereRadius;
  float cylinderRadius;
  float cylinderSpacing;
  float cylinderHeight;
  float bath_concentration;
  float sheath_concentration;
  float core_concentration;
  float D_bath;
  float D_sheath;
  float D_core;
  float gradient_strength;
  float snr;
  float dpfgAngle;
  float EulerAngle1;
  float EulerAngle2;
  float EulerAngle3;
  float ncEulerAngle1;
  float ncEulerAngle2;
  float ncEulerAngle3;
  float P_scale;
  float ncCurrent;
  int periodic;
  int parallelPlates;
  int nc_Nx;
  int nc_Ny;
  int nc_Nz;
  int repeat_number;
  int repeat_time;
  int flip_angle;
  int use_relax;
  int echo_time;
  int echo_spacing;
  int echo_number;
  int tensorFileDPFG;
  bc_t BCx0;
  bc_t BCx1;
  bc_t BCy0;
  bc_t BCy1;
  bc_t BCz0;
  bc_t BCz1;
};

/* Parameter input stuff */
struct conversion_failure { };

template <typename T>
T from_string (const std::string & s) {
  T result;
  std::istringstream stream (s);
  if (stream >> result) return result;
  throw conversion_failure ();
}

static struct {
  int t;
  int tp;
} B0 = {0,0};

void parse_input(input_data *d, char *file_name) {
  ConfigFile config(file_name);
  d->debugMCell = config.read("debug mcell", 0);
  d->Delta = config.read("Delta", 100000);
  d->delta = config.read("delta", 10000);
  d->gradient_strength = config.read("gradient strength", 4.0);
  d->bath_concentration = config.read("bath concentration", 0.95);
  d->sheath_concentration = config.read("sheath concentration", 0.5);
  d->core_concentration = config.read("core concentration", 0.88);
  d->D_bath = config.read("bath diffusion", 2.0e-5);
  d->D_sheath = config.read("sheath diffusion", 3.0e-7);
  d->D_core = config.read("core diffusion", 7.5e-6);
  d->ramp = config.read("ramp", 500);
  d->pulse = config.read("pulse", 1);
  d->mixing_time = config.read("mixing time", 0);
  d->dpfgAngle = config.read("DPFG angle", 0.0);
  d->tensorFileDPFG = config.read("tensor file DPFG", 0);
  d->EulerAngle1 = config.read("Euler rotation angle 1", 0.0);
  d->EulerAngle2 = config.read("Euler rotation angle 2", 0.0);
  d->EulerAngle3 = config.read("Euler rotation angle 3", 0.0);
  d->ncEulerAngle1 = config.read("ncCurrent Euler rotation angle 1", 0.0);
  d->ncEulerAngle2 = config.read("ncCurrent Euler rotation angle 2", 0.0);
  d->ncEulerAngle3 = config.read("ncCurrent Euler rotation angle 3", 0.0);
  d->periodic = config.read("periodic", 0);
  d->parallelPlates = config.read("parallel plates", 0);
  d->snr = config.read("snr", 100);
  d->voxel_size = config.read("voxel size", 50.0);
  d->xDim = config.read("x dim", 0.0);
  d->yDim = config.read("y dim", 0.0);
  d->zDim = config.read("z dim", 0.0);
  d->xScale = config.read("x scale", 1.0);
  d->yScale = config.read("y scale", 1.0);
  d->zScale = config.read("z scale", 1.0);
  d->volume_size = config.read("volume size", 50.0);
  d->volumeHeight = config.read("volume height", 0.0);
  d->coreRadius = config.read("core radius", 0.0);
  d->sheathRadius = config.read("sheath radius", 0.0);
  d->sphereRadius = config.read("sphere radius", 0.0);
  d->cylinderRadius = config.read("cylinder radius", 0.0);
  d->cylinderSpacing = config.read("cylinder spacing", 0.0);
  d->cylinderHeight = config.read("cylinder height", 0.0);
  d->nsteps = config.read("nsteps", 1);
  d->step_size = config.read("step size", 1.);
  d->ntess = config.read("ntess", 0);
  d->gradientDirections = config.read("gradient directions", 1);
  d->number_of_particles = config.read("nparts", 10000);
  d->randomSeed = config.read("random seed", 1);
  d->threeRegion = config.read("three region", 0);
  d->permeable = config.read("permeable", 1);
  d->P_scale = config.read("permeability scale", 1.0);
  d->ncCurrent = config.read("ncCurrent strength", 0.0);
  d->cylinderArray = config.read("cylinder array", 0);
  d->structureArray = config.read("structure array", 0);
  d->structurePacking = (structure_t)config.read("structure packing", 0);
  d->nc_Nx = config.read("nc Nx", 256);
  d->nc_Ny = config.read("nc Ny", 256);
  d->nc_Nz = config.read("nc Nz", 256);
  d->repeat_number = config.read("repeat number", 1); 
  d->repeat_time = config.read("TR", 0);
  d->flip_angle = config.read("flip angle", 90);
  d->use_relax = config.read("relaxation", 0);
  d->echo_time = config.read("TE", 0);
  d->echo_spacing = config.read("ESP", 0);
  d->echo_number = config.read("NE", 1);
  if ( !config.readInto<std::string>(d->signal_filename, "signal file",
                                     "signal_file.dat") ) {
    cout << "Using default value for signal file" << endl;
  }
  config.readInto<std::string>(d->complex_signal_filename, "complex signal file", "");
  config.readInto<std::string>(d->mcell_file, "mcell file", "");
  config.readInto<std::string>(d->log_file, "log file", "");
  config.readInto<std::string>(d->tensor_file, "tensor file", "");
  config.readInto<std::string>(d->structureFile, "structure file", "");
  config.readInto<std::string>(d->sheathFile, "sheath file", "");
  if ((d->structurePacking == NONE) && (d->cylinderArray == 1)) {
    d->structurePacking = HEX2D;
  }
  config.readInto<std::string>(d->ncStrengthFile, "ncCurrent strength file", "");
  config.readInto<std::string>(d->ncOutputFile, "ncCurrent output file", "");
  config.readInto<std::string>(d->ncObjectName, "ncCurrent object name", "");
  config.readInto<std::string>(d->ncJinFile, "ncCurrent J input file", "");
  config.readInto<std::string>(d->ncJoutFile, "ncCurrent J output file", "");
  config.readInto<std::string>(d->dpfgAngleFile, "DPFG angle file", "");
  d->BCx0 = (bc_t)config.read("BC x0", 0);
  d->BCx1 = (bc_t)config.read("BC x1", 0);
  d->BCy0 = (bc_t)config.read("BC y0", 0);
  d->BCy1 = (bc_t)config.read("BC y1", 0);
  d->BCz0 = (bc_t)config.read("BC z0", 0);
  d->BCz1 = (bc_t)config.read("BC z1", 0);
}

/* some useful debugging stuff */
void diffusionStatistics(mcellSim sim, int currentTime) {
  int rank = 0;

  /* write out a bunch of statistics */
  int i = 0;
  int numberOut = 0;
  double max_diffusion = 0.0;
  double max_radial_diffusion = 0.0;
  double max_axial_diffusion = 0.0;

  double max_ext_diffusion = 0.0;
  double max_ext_radial_diffusion = 0.0;
  double max_ext_axial_diffusion = 0.0;

  double max_x = 0.0;
  double max_y = 0.0;
  double max_z = 0.0;
  double max_ext_x = 0.0;
  double max_ext_y = 0.0;
  double max_ext_z = 0.0;
  int max_int_steps = 0;
  int max_ext_int_steps = 0;
  int core_molecules = 0;
  int sheath_molecules = 0;
  int bath_molecules = 0;
  double mean_abs_z_escaped = 0;

  int hist_size = 100;
  int unbinned_x = 0, unbinned_y = 0, unbinned_z = 0;
  gtb_histogram *hx = gtb_histogram_alloc(hist_size);
  gtb_histogram *hy = gtb_histogram_alloc(hist_size);
  gtb_histogram *hz = gtb_histogram_alloc(hist_size);
  gtb_histogram_set_uniform_ranges(hx, -12.0*sqrt(7.5e-4 * currentTime), 12.0*sqrt(7.5e-4 * currentTime));
  gtb_histogram_set_uniform_ranges(hy, -12.0*sqrt(7.5e-4 * currentTime), 12.0*sqrt(7.5e-4 * currentTime));
  gtb_histogram_set_uniform_ranges(hz, -12.0*sqrt(7.5e-4 * currentTime), 12.0*sqrt(7.5e-4 * currentTime));
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            double dir[3];
            dir[0] = 0.01*(mp->pos.x - mp->orig_pos.x);
            dir[1] = 0.01*(mp->pos.y - mp->orig_pos.y);
            dir[2] = 0.01*(mp->pos.z - mp->orig_pos.z);
            unbinned_x += gtb_histogram_increment(hx, dir[0]);
            unbinned_y += gtb_histogram_increment(hy, dir[1]);
            unbinned_z += gtb_histogram_increment(hz, dir[2]);
            float blx = dir[0] + sim.xDim;
            float bly = dir[1] + sim.yDim;
            float brx = dir[0] - sim.xDim;
            float bry = dir[1] + sim.yDim;
            float tlx = dir[0] + sim.xDim;
            float tly = dir[1] - sim.yDim;
            float trx = dir[0] - sim.xDim;
            float tRy = dir[1] - sim.yDim;
            float r2 = sim.cylinderRadius*sim.cylinderRadius;
            if (DOT_VEC3(dir, dir) > max_diffusion) {
              max_diffusion = DOT_VEC3(dir, dir);
            }
#ifdef DIFF_SIM_DEBUG
            /* for now, numberOut is really number in */
            /*    if (dir[0]*dir[0] + dir[1]*dir[1] < r2
                  || blx*blx + bly*bly < r2 || brx*brx + bry*bry < r2
                  || tlx*tlx + tly*tly < r2 || trx*trx + tRy*tRy < r2) */
            /* need to scale [xyz]Dim by 100 to get same units */
            if (fabs(mp->pos.x) > sim.xDim*100.0/2.0 || fabs(mp->pos.y) > sim.yDim*100.0/2.0 || fabs(mp->pos.z) > sim.zDim*100.0/2.0) {
              numberOut++;
              mean_abs_z_escaped += fabs(mp->pos.z);
            }
#endif /* DIFF_SIM_DEBUG */
            if (dir[0]*dir[0] + dir[1]*dir[1] > max_radial_diffusion) {
              max_radial_diffusion = dir[0]*dir[0] + dir[1]*dir[1];
            }
            if (fabs(dir[2]) > max_axial_diffusion) {
              max_axial_diffusion = fabs(dir[2]);
            }
            if (fabs(dir[0]) > max_x) {
              max_x = fabs(dir[0]);
            }
            if (fabs(dir[1]) > max_y) {
              max_y = fabs(dir[1]);
            }
            if (fabs(dir[2]) > max_z) {
              max_z = fabs(dir[2]);
            }
          }
        }
      }
    }
  }


  if (rank == 0) {
    /*char hfilenameX[256];
    char hfilenameY[256];
    char hfilenameZ[256];
    snprintf(hfilenameX, 256, "hx-%d.dat", currentTime);
    snprintf(hfilenameY, 256, "hy-%d.dat", currentTime);
    snprintf(hfilenameZ, 256, "hz-%d.dat", currentTime);
    ofstream hfileX(hfilenameX, ios::out);
    ofstream hfileY(hfilenameY, ios::out);
    ofstream hfileZ(hfilenameZ, ios::out);
    for (int i = 0; i < hx->size; i++) {
      hfileX << gtb_histogram_get_bin_mean(hx, i) << " "
             << gtb_histogram_get_value(hx, i) << endl;
    }
    for (int i = 0; i < hy->size; i++) {
      hfileY << gtb_histogram_get_bin_mean(hy, i) << " "
             << gtb_histogram_get_value(hy, i) << endl;
    }
    for (int i = 0; i < hz->size; i++) {
      hfileZ << gtb_histogram_get_bin_mean(hz, i) << " "
             << gtb_histogram_get_value(hz, i) << endl;
    }
    hfileX.close();
    hfileY.close();
    hfileZ.close();*/
    printf("Maximum internal total diffusion: %f microns\n", sqrt(max_diffusion));
    printf("Maximum internal radial diffusion: %f microns\n", sqrt(max_radial_diffusion));
    printf("Maximum internal axial diffusion: %f microns\n", max_axial_diffusion);
    printf("max x, y, z diffusion: [%f, %f, %f]\n", max_x, max_y, max_z);
#ifdef DIFF_SIM_DEBUG
    if (numberOut > 0) {
      printf("abs z location of escapees: %f\n", mean_abs_z_escaped/numberOut);
    }
    printf("%d molecules are outside the cylinder radius\n", numberOut);
#endif /* DIFF_SIM_DEBUG */
  }
  gtb_histogram_free(hx);
  gtb_histogram_free(hy);
  gtb_histogram_free(hz);
}
void diffusionStatistics(mcellSim sim, double &max_x, double &max_y,
                         double &max_z) {
  int i = 0;
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  max_x = max_y = max_z = 0;
  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            double dir[3];
            dir[0] = 0.01*(mp->pos.x - mp->orig_pos.x);
            dir[1] = 0.01*(mp->pos.y - mp->orig_pos.y);
            dir[2] = 0.01*(mp->pos.z - mp->orig_pos.z);
            if (fabs(dir[0]) > max_x) {
              max_x = fabs(dir[0]);
            }
            if (fabs(dir[1]) > max_y) {
              max_y = fabs(dir[1]);
            }
            if (fabs(dir[2]) > max_z) {
              max_z = fabs(dir[2]);
            }
          }
        }
      }
    }
  }
}

/* flip sign of phase to implement 180 RF pulse */
void rf_180(int nDirs) {
  int i, j;
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            for (j = 0; j < nDirs; j++) {
              mp->phase_shift[j] *= -1.0;
            }
          }
        }
      }
    }
  }
}

void rf_pulse(double alpha, int nDirs, int aboutaxis, int rand) {
  //cout << "RF pulse - flip angle " << alpha << "\n";
  double sa = sin(M_PI*alpha/180.);
  double ca = cos(M_PI*alpha/180.);
  double Mxm,Mym,Mzm;
  int i, j;
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            if (rand) {
              alpha = 360.*wrap_rng();
              sa = sin(M_PI*alpha/180.);
              ca = cos(M_PI*alpha/180.);
            }
            for (j = 0; j < nDirs; j++) {
              int dt = B0.t-B0.tp;
              if (dt>0) {
                mp->Mx[j] *= exp(-(double)dt/mp->properties->relax_t2);
                mp->My[j] *= exp(-(double)dt/mp->properties->relax_t2);
                mp->Mz[j] = mp->Mz[j]*exp(-(double)dt/mp->properties->relax_t1) + (1 - exp(-(double)dt/mp->properties->relax_t1));
              }
              Mxm = mp->Mx[j];
              Mym = mp->My[j];
              Mzm = mp->Mz[j];
              if (aboutaxis == 0) {
                mp->Mx[j] = Mxm;
                mp->My[j] = Mym*ca - Mzm*sa;
                mp->Mz[j] = Mym*sa + Mzm*ca;
              } else if (aboutaxis == 1) {
                mp->Mx[j] = Mxm*ca + Mzm*sa;
                mp->My[j] = Mym;
                mp->Mz[j] = -Mxm*sa + Mzm*ca;
              } else if (aboutaxis == 2) {
                mp->Mx[j] = Mxm*ca - Mym*sa;
                mp->My[j] = Mxm*sa + Mym*ca;
                mp->Mz[j] = Mzm;
              }
              mp->phase_shift[j] = atan2(mp->My[j],mp->Mx[j]);
            }
          }
        }
      }
    }
  }
  B0.tp = B0.t;
}

/* create the mcell file template */
char *mcellTemplate(string structure, string sheath, mcellSim sim, int nProcessors) {
  char *mcellFileName = mktemp(strdup("mcell.XXXXXX"));
  ofstream mcellFile(mcellFileName, ios::out);
  int threeRegion = sim.threeRegion;
  int permeable = sim.permeable;

  mcellFile << "dt = " << sim.timestep * 1.0e-6 << endl
            << "sgd = 10000  /* surface grid density */" << endl
            << "ITERATIONS = " << sim.endTime/sim.timestep << endl
            << "TIME_STEP = dt" << endl
            << "TIME_STEP_MAX = dt" << endl
            << "SURFACE_GRID_DENSITY = sgd" << endl
            << "Dwc = " << sim.D_core << endl
            << "P_scale = " << sim.P_scale << endl
            << "Na = 6.0221417930e23" << endl
            << "concentration = 100*1e15/Na" << endl
            << "conc_c = " << sim.core_concentration << "*concentration" << endl 
            << "MICROSCOPIC_REVERSIBILITY = TRUE" << endl;
  if (threeRegion) {
    mcellFile << "Dwb = " << sim.D_bath << endl
              << "Dws = " << sim.D_sheath << endl
              << "conc_b = " << sim.bath_concentration << "*concentration" << endl
              << "conc_s = " << sim.sheath_concentration << "*concentration" << endl<< endl;
    if (sim.core_concentration < EPS_C) {
      mcellFile << "p_cs = 1.0" << endl;
    }
    else {
      mcellFile << "p_cs = SQRT(Dws/Dwc)*conc_s/conc_c" << endl;
    }
    if (sim.bath_concentration < EPS_C) {
      mcellFile << "p_bs = 1.0" << endl << endl;
    }
    else {
      mcellFile << "p_bs = SQRT(Dws/Dwb)*conc_s/conc_b" << endl << endl;
    }
    mcellFile << "p_cs_fix = P_Scale*p_cs" << endl
              << "p_bs_fix = P_Scale*p_bs" << endl << endl;
    if (permeable) {
      mcellFile << "rate_core_sheath = p_cs_fix*(Na/sgd)*SQRT(1e8*Dwc/(PI*dt))/1e15" << endl
                << "rate_sheath_core = P_scale*(Na/sgd)*SQRT(1e8*Dws/(PI*dt))/1e15" << endl
                << "rate_sheath_bath = rate_sheath_core" << endl
                << "rate_bath_sheath = p_bs_fix*(Na/sgd)*SQRT(1e8*Dwb/(PI*dt))/1e15" << endl << endl;
    }
    mcellFile << "DEFINE_SURFACE_CLASSES {" << endl
              << "  interface_core_sheath {" << endl
              << "    REFLECTIVE = ALL_MOLECULES" << endl
              << "  }" << endl
              << "  interface_sheath_bath {" << endl
              << "    REFLECTIVE = ALL_MOLECULES" << endl
              << "  }" << endl;
    if (sim.periodic) {
      mcellFile << "  periodic_boundary {" << endl
                << "    PERIODIC = ALL_MOLECULES" << endl
                << "  }" << endl;
    }
    mcellFile << "}" << endl;
  }

  if (sim.xDim == 0.0) {
    mcellFile << "PARTITION_X = [[" << -1.01*sim.volume_size/2.0
              << " TO " << 1.01*sim.volume_size/2.0
              << " STEP 0.123456789]]" << endl;
  }
  else {
    mcellFile << "PARTITION_X = [[" << -1.01*sim.xDim/2.0
              << " TO " << 1.01*sim.xDim/2.0
              << " STEP 0.123456789]]" << endl;
  }
  if (sim.yDim == 0.0) {
    mcellFile << "PARTITION_Y = [[" << -1.01*sim.volume_size/2.0
              << " TO " << 1.01*sim.volume_size/2.0
              << " STEP 0.123456789" << "]]" << endl;
  }
  else {
    mcellFile << "PARTITION_Y = [[" << -1.01*sim.yDim/2.0
              << " TO " << 1.01*sim.yDim/2.0
              << " STEP 0.123456789" << "]]" << endl;
  }
  if (sim.zDim == 0.0) {
    mcellFile << "PARTITION_Z = [[" << -1.01*sim.volume_size/2.0
              << " TO " << 1.01*sim.volume_size/2.0
              << " STEP 0.123456789" << "]]" << endl;
  }
  else {
    mcellFile << "PARTITION_Z = [[" << -1.01*sim.zDim/2.0
              << " TO " << 1.01*sim.zDim/2.0
              << " STEP 0.123456789" << "]]" << endl;
  }
  if (threeRegion) {
    mcellFile << "DEFINE_MOLECULES {" << endl
              << "  water_core {DIFFUSION_CONSTANT_3D = Dwc}" << endl
              << "  water_sheath {DIFFUSION_CONSTANT_3D = Dws}" << endl
              << "  water_bath {DIFFUSION_CONSTANT_3D = Dwb}" << endl
              << "}" << endl << endl;
    if (permeable) {
      mcellFile << "DEFINE_REACTIONS {" << endl
                << "  water_core, @ interface_core_sheath' -> water_sheath' [rate_core_sheath] : core_sheath" << endl
                << "  water_sheath' @ interface_core_sheath' -> water_core, [rate_sheath_core] : sheath_core" << endl
                << "  water_sheath,, @ interface_sheath_bath'' -> water_bath'' [rate_sheath_bath] : sheath_bath" << endl
                << "  water_bath'' @ interface_sheath_bath'' -> water_sheath,, [rate_bath_sheath] : bath_sheath" << endl
                << "}" << endl << endl;
    }
  }
  else {
    mcellFile << "DEFINE_MOLECULES {" << endl
      /* water at room temperature */
              << " water {DIFFUSION_CONSTANT_3D = Dwc}" << endl
#if 1
              << " water_internal {DIFFUSION_CONSTANT_3D = Dwc}" << endl
              << " water_external {DIFFUSION_CONSTANT_3D = Dwc}" << endl
#endif
              << "}" << endl << endl;
  }

  if (sim.periodic && !threeRegion) {
    mcellFile << "DEFINE_SURFACE_CLASSES {" << endl
              << "  membrane {" << endl
              << "    REFLECTIVE = ALL_MOLECULES" << endl
              << "  }" << endl
              << "  periodic_boundary {" << endl
              << "    PERIODIC = ALL_MOLECULES" << endl
              << "  }" << endl << "}" << endl;
  }

  if (structure != "") {
    mcellFile << "INCLUDE_FILE = \"" << structure << "\"" << endl << endl;
  }
  if (sheath != "") {
    mcellFile << "INCLUDE_FILE = \"" << sheath << "\"" << endl << endl;
  }

  if (permeable) {
    mcellFile << "permeability = P_scale*(Na/sgd)*SQRT(1e8*Dwc/(PI*dt))/1e15" << endl;
    mcellFile << "DEFINE_REACTIONS {" << endl
              << "  water; @ membrane; -> water; [permeability] : water_water" << endl
#if 1
              << "  water_internal, @ membrane; -> water_external' [permeability] : internal_external" << endl
              << "  water_external, @ membrane; -> water_internal' [permeability] : external_internal" << endl
#endif
              << "}" << endl << endl;
  }

  float boxXMin = (sim.xDim == 0.0) ? -sim.volume_size/2.0 : -sim.xDim/2.0;
  float boxYMin = (sim.yDim == 0.0) ? -sim.volume_size/2.0 : -sim.yDim/2.0;
  float boxZMin = (sim.zDim == 0.0) ? -sim.volume_size/2.0 : -sim.zDim/2.0;
  float boxXMax = (sim.xDim == 0.0) ? sim.volume_size/2.0 : sim.xDim/2.0;
  float boxYMax = (sim.yDim == 0.0) ? sim.volume_size/2.0 : sim.yDim/2.0;
  float boxZMax = (sim.zDim == 0.0) ? sim.volume_size/2.0 : sim.zDim/2.0;

  mcellFile << "boundingBox BOX {" << endl
            << "  CORNERS = [" 
            << boxXMin << " , " << boxYMin << ", " << boxZMin << "], ["
            << boxXMax << " , " << boxYMax << ", " << boxZMax << "]" << endl;

  if (sim.periodic) {
    if (sim.parallelPlates) {
      mcellFile << "  DEFINE_SURFACE_REGIONS {" << endl
                << "    interface {" << endl
                << "      ELEMENT_LIST = [TOP, BOTTOM, FRONT, BACK]" << endl
                << "      SURFACE_CLASS = periodic_boundary" << endl
                << "    }" << endl << "  }" << endl;
    }
    else {
      mcellFile << "  DEFINE_SURFACE_REGIONS {" << endl
                << "    interface {" << endl
                << "      ELEMENT_LIST = [ALL_ELEMENTS]" << endl
                << "      SURFACE_CLASS = periodic_boundary" << endl
                << "    }" << endl << "  }" << endl;
    }
  }
  mcellFile << "}" << endl << endl;

  mcellFile << "INSTANTIATE diffsim_domain OBJECT {" << endl;

  mcellFile << "  boundingBoxObject OBJECT boundingBox {}" << endl;

  float radialScale = 1.0;
  float lengthScale = 1.0;
  if (sim.structurePacking == HEX2D) {
    if (sim.coreRadius > 0.0 || sim.cylinderRadius > 0.0) {
      radialScale = max(sim.coreRadius, sim.cylinderRadius);
    }
    if (sim.cylinderHeight > 0.0) {
      lengthScale = sim.cylinderHeight;
    }
  }
  if (structure != "") {
    if ((sim.structurePacking == HEX2D) && sim.xDim != 0.0 && sim.yDim != 0.0) {
      if (threeRegion) {
        /* Cylinder array with core-sheath-bath regions */

        /* center */
        mcellFile << "  centerCoreObject OBJECT coreStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.coreRadius << ", " << 2.0*sim.coreRadius << ", " << lengthScale << "]" << endl
                  << "  }" << endl
                  << "  centerSheathObject OBJECT sheathStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.sheathRadius << ", " << 2.0*sim.sheathRadius << ", " << 1.1*lengthScale << "]" << endl
                  << "  }" << endl

          /* top left */
                  << "  topLeftCoreObject OBJECT coreStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.coreRadius << ", " << 2.0*sim.coreRadius << ", " << lengthScale << "]" << endl
                  << "    TRANSLATE = [" << -sim.xDim/2.0 << ", " 
                  << sim.yDim/2.0 << ", 0]" << endl
                  << "  }" << endl
                  << "  topLeftSheathObject OBJECT sheathStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.sheathRadius << ", " << 2.0*sim.sheathRadius << ", " << 1.1*lengthScale << "]" << endl
                  << "    TRANSLATE = [" << -sim.xDim/2.0 << ", " 
                  << sim.yDim/2.0 << ", 0]" << endl
                  << "  }" << endl

          /* top right */
                  << "  topRightCoreObject OBJECT coreStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.coreRadius << ", " << 2.0*sim.coreRadius << ", " << lengthScale << "]" << endl
                  << "    TRANSLATE = [" << sim.xDim/2.0 << ", " 
                  << sim.yDim/2.0 << ", 0]" << endl
                  << "  }" << endl
                  << "  topRightSheathObject OBJECT sheathStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.sheathRadius << ", " << 2.0*sim.sheathRadius << ", " << 1.1*lengthScale << "]" << endl
                  << "    TRANSLATE = [" << sim.xDim/2.0 << ", " 
                  << sim.yDim/2.0 << ", 0]" << endl
                  << "  }" << endl

          /* bottom left */
                  << "  bottomLeftCoreObject OBJECT coreStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.coreRadius << ", " << 2.0*sim.coreRadius << ", " << lengthScale << "]" << endl
                  << "    TRANSLATE = [" << -sim.xDim/2.0 << ", " 
                  << -sim.yDim/2.0 << ", 0]" << endl
                  << "  }" << endl
                  << "  bottomLeftSheathObject OBJECT sheathStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.sheathRadius << ", " << 2.0*sim.sheathRadius << ", " << 1.1*lengthScale << "]" << endl
                  << "    TRANSLATE = [" << -sim.xDim/2.0 << ", " 
                  << -sim.yDim/2.0 << ", 0]" << endl
                  << "  }" << endl

          /* bottom right */
                  << "  bottomRightCoreObject OBJECT coreStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.coreRadius << ", " << 2.0*sim.coreRadius << ", " << lengthScale << "]" << endl
                  << "    TRANSLATE = [" << sim.xDim/2.0 << ", " 
                  << -sim.yDim/2.0 << ", 0]" << endl
                  << "  }" << endl
                  << "  bottomRightSheathObject OBJECT sheathStructure {" << endl
                  << "    SCALE = [" << 2.0*sim.sheathRadius << ", " << 2.0*sim.sheathRadius << ", " << 1.1*lengthScale << "]" << endl
                  << "    TRANSLATE = [" << sim.xDim/2.0 << ", " 
                  << -sim.yDim/2.0 << ", 0]" << endl
                  << "  }" << endl;

        /* release sites */
        mcellFile << "  coreReleaseSite RELEASE_SITE {" << endl
                  << "    SHAPE = diffsim_domain.boundingBoxObject[ALL]*(diffsim_domain.topLeftCoreObject[ALL] + diffsim_domain.topRightCoreObject[ALL] + diffsim_domain.bottomLeftCoreObject[ALL] + diffsim_domain.bottomRightCoreObject[ALL] + diffsim_domain.centerCoreObject[ALL])"
                  << "    MOLECULE = water_core" << endl
                  << "    CONCENTRATION = conc_c" << endl
                  << "  }" << endl;
        mcellFile << "  sheathReleaseSite RELEASE_SITE {" << endl
                  << "    SHAPE = diffsim_domain.boundingBoxObject[ALL]*(diffsim_domain.topLeftSheathObject[ALL] + diffsim_domain.topRightSheathObject[ALL] + diffsim_domain.bottomLeftSheathObject[ALL] + diffsim_domain.bottomRightSheathObject[ALL] + diffsim_domain.centerSheathObject[ALL] - diffsim_domain.topLeftCoreObject[ALL] - diffsim_domain.topRightCoreObject[ALL] - diffsim_domain.bottomLeftCoreObject[ALL] - diffsim_domain.bottomRightCoreObject[ALL] - diffsim_domain.centerCoreObject[ALL])"
                  << "    MOLECULE = water_sheath" << endl
                  << "    CONCENTRATION = conc_s" << endl
                  << "  }" << endl;
        mcellFile << "  bathReleaseSite RELEASE_SITE {" << endl
                  << "    SHAPE = diffsim_domain.boundingBoxObject[ALL] - (diffsim_domain.topLeftSheathObject[ALL] + diffsim_domain.topRightSheathObject[ALL] + diffsim_domain.bottomLeftSheathObject[ALL] + diffsim_domain.bottomRightSheathObject[ALL] + diffsim_domain.centerSheathObject[ALL])"
                  << "    MOLECULE = water_bath" << endl
                  << "    CONCENTRATION = conc_b" << endl
                  << "  }" << endl;
      }
      else {
        /* 2d hex packed cylinder array with uniform diffusion */
        mcellFile << "  topLeftObject OBJECT structure { TRANSLATE = [" << -sim.xDim/2.0 << ", " << sim.yDim/2.0 << ", 0] }" << endl;
        mcellFile << "  topRightObject OBJECT structure { TRANSLATE = [" << sim.xDim/2.0 << ", " << sim.yDim/2.0 << ", 0] }" << endl;
        mcellFile << "  bottomLeftObject OBJECT structure { TRANSLATE = [" << -sim.xDim/2.0 << ", " << -sim.yDim/2.0 << ", 0] }" << endl;
        mcellFile << "  bottomRightObject OBJECT structure { TRANSLATE = [" << sim.xDim/2.0 << ", " << -sim.yDim/2.0 << ", 0] }" << endl;
        mcellFile << "  centerObject OBJECT structure {}" << endl;
        int Nparts = (int)(sim.num_particles*sim.cylinderRadius*sim.cylinderRadius*
                           M_PI*2.0/sim.xDim/sim.yDim);
#if 1
        mcellFile << "  internalReleaseSite RELEASE_SITE {" << endl
                  << "    SHAPE = diffsim_domain.boundingBoxObject[ALL]*(diffsim_domain.topLeftObject[ALL] + diffsim_domain.topRightObject[ALL] + diffsim_domain.bottomLeftObject[ALL] + diffsim_domain.bottomRightObject[ALL] + diffsim_domain.centerObject[ALL])"
                  << "    MOLECULE = water_internal" << endl
                  << "    NUMBER_TO_RELEASE = " << Nparts/nProcessors << endl
                  << "  }" << endl;
        mcellFile << "  externalReleaseSite RELEASE_SITE {" << endl
                  << "    SHAPE = diffsim_domain.boundingBoxObject[ALL]-(diffsim_domain.topLeftObject[ALL] + diffsim_domain.topRightObject[ALL] + diffsim_domain.bottomLeftObject[ALL] + diffsim_domain.bottomRightObject[ALL] + diffsim_domain.centerObject[ALL])"
                  << "    MOLECULE = water_external" << endl
                  << "    NUMBER_TO_RELEASE = " << (sim.num_particles-Nparts)/nProcessors << endl
                  << "  }" << endl;
#else
        mcellFile << "  waterReleaseSite RELEASE_SITE {" << endl
                  << "    SHAPE = diffsim_domain.boundingBoxObject[ALL]" << endl
                  << "    MOLECULE = water" << endl
                  << "    NUMBER_TO_RELEASE = " << sim.num_particles/nProcessors << endl
                  << "  }" << endl;
#endif
      }
    } else if ( (sim.structurePacking == CCP3D) && 
                (sim.xDim != 0.0) && (sim.yDim != 0.0) && (sim.zDim != 0.0) ) {
        float xs = sim.xScale;
        float ys = sim.yScale;
        float zs = sim.zScale;
        /* 3D cubic close packed array with uniform diffusion */
        mcellFile << "  tObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << 0.0 << ", " << 0.0 << ", " 
                  << -sim.zDim/2.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  bObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << 0.0 << ", " << 0.0 << ", " 
                  << sim.zDim/2.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  lObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << 0.0 << ", " << -sim.yDim/2.0 
                  << ", " << 0.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  rObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << 0.0 << ", " << sim.yDim/2.0 
                  << ", " << 0.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  pObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << -sim.xDim/2.0 << ", " << -0.0 
                  << ", " << 0.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  aObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << sim.xDim/2.0 << ", " << -0.0 
                  << ", " << 0.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  tlpObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << -sim.xDim/2.0 << ", " 
                  << -sim.yDim/2.0 << ", " << -sim.zDim/2.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  tlaObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" <<  sim.xDim/2.0 << ", " 
                  << -sim.yDim/2.0 << ", " << -sim.zDim/2.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  trpObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << -sim.xDim/2.0 << ", " 
                  <<  sim.yDim/2.0 << ", " << -sim.zDim/2.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  traObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" <<  sim.xDim/2.0 << ", " 
                  <<  sim.yDim/2.0 << ", " << -sim.zDim/2.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  blpObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << -sim.xDim/2.0 << ", " 
                  << -sim.yDim/2.0 << ", " <<  sim.zDim/2.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  blaObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" <<  sim.xDim/2.0 << ", " 
                  << -sim.yDim/2.0 << ", " <<  sim.zDim/2.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        mcellFile << "  brpObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" << -sim.xDim/2.0 << ", " 
                  <<  sim.yDim/2.0 << ", " <<  sim.zDim/2.0 << "]"  << endl
                  << "  }" << endl;
        mcellFile << "  braObject OBJECT structure {" << endl 
                  << "    TRANSLATE = [" <<  sim.xDim/2.0 << ", " 
                  <<  sim.yDim/2.0 << ", " <<  sim.zDim/2.0 << "]"  << endl
                  << "    SCALE = [" << xs << ", " << ys << ", " << zs << "]" 
                  << "  }" << endl;
        int Nparts = (int)(sim.num_particles*M_PI*16.0/3.0*sim.sphereRadius*
                           sim.sphereRadius*sim.sphereRadius/
                           sim.xDim/sim.yDim/sim.zDim);
        mcellFile << "  internalReleaseSite RELEASE_SITE {" << endl
                  << "    SHAPE = diffsim_domain.boundingBoxObject[ALL]*(diffsim_domain.tObject[ALL] + diffsim_domain.bObject[ALL] + diffsim_domain.lObject[ALL] + diffsim_domain.rObject[ALL] + diffsim_domain.aObject[ALL]+diffsim_domain.pObject[ALL] + diffsim_domain.tlpObject[ALL] + diffsim_domain.tlaObject[ALL] + diffsim_domain.trpObject[ALL] + diffsim_domain.traObject[ALL] + diffsim_domain.blpObject[ALL] + diffsim_domain.blaObject[ALL] + diffsim_domain.brpObject[ALL] + diffsim_domain.braObject[ALL])"
                  << "    MOLECULE = water_internal" << endl
                  << "    NUMBER_TO_RELEASE = " << Nparts/nProcessors << endl
                  << "  }" << endl;
        mcellFile << "  externalReleaseSite RELEASE_SITE {" << endl
                  << "    SHAPE = diffsim_domain.boundingBoxObject[ALL]-(diffsim_domain.tObject[ALL] + diffsim_domain.bObject[ALL] + diffsim_domain.lObject[ALL] + diffsim_domain.rObject[ALL] + diffsim_domain.aObject[ALL]+diffsim_domain.pObject[ALL] + diffsim_domain.tlpObject[ALL] + diffsim_domain.tlaObject[ALL] + diffsim_domain.trpObject[ALL] + diffsim_domain.traObject[ALL] + diffsim_domain.blpObject[ALL] + diffsim_domain.blaObject[ALL] + diffsim_domain.brpObject[ALL] + diffsim_domain.braObject[ALL])"
                  << "    MOLECULE = water_external" << endl
                  << "    NUMBER_TO_RELEASE = " << (sim.num_particles-Nparts)/nProcessors << endl
                  << "  }" << endl;
    }
    else {
      mcellFile << "  structureObject OBJECT structure {" << endl
                << "    SCALE = [" << radialScale << ", " << radialScale << ", " << lengthScale << "]" << endl
                << "  }" << endl;

      mcellFile << "  waterReleaseSite RELEASE_SITE {" << endl
                << "    SHAPE = diffsim_domain.structureObject[ALL] * diffsim_domain.boundingBoxObject[ALL]" << endl
                << "    MOLECULE = water" << endl
                << "    NUMBER_TO_RELEASE = " << sim.num_particles/nProcessors << endl
                << "  }" << endl;
    }
  }
  else {
    mcellFile << "  waterReleaseSite RELEASE_SITE {" << endl
              << "    SHAPE = diffsim_domain.boundingBoxObject[ALL]" << endl
              << "    MOLECULE = water" << endl
              << "    NUMBER_TO_RELEASE = " << sim.num_particles/nProcessors << endl
              << "  }" << endl;
  }

  mcellFile << "}" << endl;

  if (threeRegion) {
    mcellFile << "REACTION_DATA_OUTPUT {\n\
  STEP = dt\n\
  {COUNT[water_core,WORLD]} => \"h2o_c.world.dat\"\n\
  {COUNT[water_sheath,WORLD]} => \"h2o_s.world.dat\"\n\
  {COUNT[water_bath,WORLD]} => \"h2o_b.world.dat\"\n\
  {COUNT[water_core,diffsim_domain.centerCoreObject,ALL_HITS]\n\
   +COUNT[water_core,diffsim_domain.topLeftCoreObject,ALL_HITS]\n\
   +COUNT[water_core,diffsim_domain.topRightCoreObject,ALL_HITS]\n\
   +COUNT[water_core,diffsim_domain.bottomLeftCoreObject,ALL_HITS]\n\
   +COUNT[water_core,diffsim_domain.bottomRightCoreObject,ALL_HITS]\n\
  } => \"h2o_c.hits.dat\"\n\
  {COUNT[water_sheath,diffsim_domain.centerCoreObject,ALL_HITS]\n\
   +COUNT[water_sheath,diffsim_domain.topLeftCoreObject,ALL_HITS]\n\
   +COUNT[water_sheath,diffsim_domain.topRightCoreObject,ALL_HITS]\n\
   +COUNT[water_sheath,diffsim_domain.bottomLeftCoreObject,ALL_HITS]\n\
   +COUNT[water_sheath,diffsim_domain.bottomRightCoreObject,ALL_HITS]\n\
  } => \"h2o_sc.hits.dat\"\n\
  {COUNT[water_sheath,diffsim_domain.centerSheathObject,ALL_HITS]\n\
   +COUNT[water_sheath,diffsim_domain.topLeftSheathObject,ALL_HITS]\n\
   +COUNT[water_sheath,diffsim_domain.topRightSheathObject,ALL_HITS]\n\
   +COUNT[water_sheath,diffsim_domain.bottomLeftSheathObject,ALL_HITS]\n\
   +COUNT[water_sheath,diffsim_domain.bottomRightSheathObject,ALL_HITS]\n\
  } => \"h2o_sb.hits.dat\"\n\
  {COUNT[water_bath,diffsim_domain.centerSheathObject,ALL_HITS]\n\
   +COUNT[water_bath,diffsim_domain.topLeftSheathObject,ALL_HITS]\n\
   +COUNT[water_bath,diffsim_domain.topRightSheathObject,ALL_HITS]\n\
   +COUNT[water_bath,diffsim_domain.bottomLeftSheathObject,ALL_HITS]\n\
   +COUNT[water_bath,diffsim_domain.bottomRightSheathObject,ALL_HITS]\n\
  } => \"h2o_b.hits.dat\"\n";
    if (permeable) {
      mcellFile << "{COUNT[core_sheath,WORLD]} => \"core_sheath.world.dat\"\n\
  {COUNT[sheath_core,WORLD]} => \"sheath_core.world.dat\"\n\
  {COUNT[sheath_bath,WORLD]} => \"sheath_bath.world.dat\"\n\
  {COUNT[bath_sheath,WORLD]} => \"bath_sheath.world.dat\"\n";
    }
    mcellFile << "}";
  }

  mcellFile.close();
  return mcellFileName;
}

/* update phase of all diffusing molecules */
void mcell_update_spins(mcellSim *s, int nDirs,
                        double (*dirs)[3], int rank) {
  int i = 0, j;
  struct ligand *curr_lig;
  double dir[3];
  double ncG[3],ncGt[3];
  double proj;
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  /* 1.0e2 (cm/m) 
   * 1.0e-6 (m/micron)
   * 26752 (gamma in rad/G/s)
   * 1.0e-6 (s/micro-s)
   = 26752e-10 rad*cm/(G*micron*micro-s) */
  /* factor is in units of rad/micron */
  double factor = s->current_gradient_strength * s->timestep * 26752e-10;
  double factorNC = s->ncCurrent * s->timestep * 4.0*M_PI*26752e-15;
  complex<double> I(0.0, 1.0);
  double rotposx, rotposy, rotposz;

  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            /* We want the actual location, not the displacement from t=0,
               but with periodic boundary conditions the actual
               location needs to be derived from the displacement. If
               the particle is currently at x and x_0 is within the
               periodic boundary conditions, then x is the correct
               location. If x_0 is outside the periodic boundary
               conditions (x_0 < x_min or x_0 > x_max) everything
               needs to be translated:
               x_0 -> x_0
               - (x_max - x_min)*floor((x_0 - x_min)/(x_max - x_min)). */
#if 0      
            dir[0] = mp->pos.x - mp->orig_pos.x;
            dir[1] = mp->pos.y - mp->orig_pos.y;
            dir[2] = mp->pos.z - mp->orig_pos.z;
            if (mp->orig_pos.x < -s->xDim/2.0 || mp->orig_pos.x > s->xDim/2.0) {
              dir[0] -= (s->xDim)*floor((mp->orig_pos.x - s->xDim/2.0)/s->xDim);
            }
            if (mp->orig_pos.y < -s->yDim/2.0 || mp->orig_pos.y > s->xDim/2.0) {
              dir[1] -= (s->yDim)*floor((mp->orig_pos.y - s->yDim/2.0)/s->yDim);
            }
            if (mp->orig_pos.z < -s->zDim/2.0 || mp->orig_pos.z > s->xDim/2.0) {
              dir[2] -= (s->zDim)*floor((mp->orig_pos.z - s->zDim/2.0)/s->zDim);
            }
#endif
          /*dir[0] = mp->pos.x + mp->offset.x;
            dir[1] = mp->pos.y + mp->offset.y;
            dir[2] = mp->pos.z + mp->offset.z;*/
            dir[0] = mp->pos.x - mp->orig_pos.x + mp->offset.x;
            dir[1] = mp->pos.y - mp->orig_pos.y + mp->offset.y;
            dir[2] = mp->pos.z - mp->orig_pos.z + mp->offset.z;

            // Rotation might need to be in opposite direction, would need to do transpose
            if (s->ncCurrent!=0) {
              rotposx = s->ncTr[0][0]*0.01*mp->pos.x + s->ncTr[0][1]*0.01*mp->pos.y + s->ncTr[0][2]*0.01*mp->pos.z;
              rotposy = s->ncTr[1][0]*0.01*mp->pos.x + s->ncTr[1][1]*0.01*mp->pos.y + s->ncTr[1][2]*0.01*mp->pos.z;
              rotposz = s->ncTr[2][0]*0.01*mp->pos.x + s->ncTr[2][1]*0.01*mp->pos.y + s->ncTr[2][2]*0.01*mp->pos.z;
              s->vp->getGradAt(rotposx, rotposy, rotposz, ncGt[0], ncGt[1], ncGt[2]);
              ncG[0] = s->ncTr[0][0]*ncGt[0] + s->ncTr[0][1]*ncGt[1] + s->ncTr[0][2]*ncGt[2];
              ncG[1] = s->ncTr[1][0]*ncGt[0] + s->ncTr[1][1]*ncGt[1] + s->ncTr[1][2]*ncGt[2];
              ncG[2] = s->ncTr[2][0]*ncGt[0] + s->ncTr[2][1]*ncGt[1] + s->ncTr[2][2]*ncGt[2];
            }
            for (j = 0; j < nDirs; j++) {
              /* MCell-3 scales units by a factor of 100. */
              proj = DOT_VEC3(dir,dirs[j]) * 0.01 * factor;
              mp->phase_shift[j] += proj;
              if (s->ncCurrent!=0) {
                proj = DOT_VEC3(dir,ncG) * 0.01 * factorNC;
                mp->phase_shift[j] += proj;
              }
              
            }
          }
        }
      }
    }
  }
}

int update_nc_current(mcellSim *s, int currentTime)
{
  int ic;
  double t0, nc;
  mcellSim &sim = *(mcellSim *)s;
  if (sim.ncCurrentTimes.size()<2) return 0;
  t0 = fmod(currentTime,(*(sim.ncCurrentTimes.end()-1)));
  ic = upper_bound(sim.ncCurrentTimes.begin(),sim.ncCurrentTimes.end(),t0) -
       sim.ncCurrentTimes.begin();
  nc = sim.ncCurrentValues[ic-1] +
       (sim.ncCurrentValues[ic] - sim.ncCurrentValues[ic-1])/
       (sim.ncCurrentTimes[ic] - sim.ncCurrentTimes[ic-1])*
       (t0 - sim.ncCurrentTimes[ic-1]);
  sim.ncCurrent = nc;
  return 0;
}

int test_point_inside(double x, double y, double z, string obj_name, mcellSim *s) {
  mcellSim &sim = *(mcellSim *)s;
  struct wall *w;
  struct wall_list *wl;
  struct subvolume *sv, *nsv, *wsv;
  struct vector3 start, end, col_hitpt, delta;
  struct waypoint *wp;
  struct subvolume *my_sv;
  int hit, col = 0;
  int allhits = 0;
  double first_time = FOREVER, col_time = 1.;
  
  start.x = 100*x;
  start.y = 100*y;
  start.z = 100*z;
  end.x = 100*sim.xDim;
  end.y = 100*sim.yDim;
  end.z = 100*sim.zDim;

  const int px = bisect_wrap(world->x_partitions,world->nx_parts,start.x);
  const int py = bisect_wrap(world->y_partitions,world->ny_parts,start.y);
  const int pz = bisect_wrap(world->z_partitions,world->nz_parts,start.z);
  const int this_sv = pz + (world->nz_parts-1)*( py + (world->ny_parts-1)*px );
  wp = &(world->waypoints[this_sv]);
  my_sv = &(world->subvol[this_sv]);
  struct region_list *rl;
  string pbox = "periodic box";
  if (wp->regions) {
    for (rl = wp->regions; rl != NULL; rl = rl->next) {
      if (pbox.compare(rl->reg->parent->last_name) > 0) {
        return 1;
      }
    }
  }
  
  for (sv = my_sv; sv != NULL; sv = next_subvol_wrap(&start,&end,sv)) {
    delta.x = wp->loc.x - start.x;
    delta.y = wp->loc.y - start.y;
    delta.z = wp->loc.z - start.z;

    double t_sv_hit = collide_sv_time_wrap(&start,&delta,sv);
    if (t_sv_hit > 1.0) t_sv_hit = 1.0;
    for (wl = sv->wall_head; wl != NULL; wl = wl->next) {
      if (wl->this_wall->num_surf_classes == 0) {
        if (obj_name == "" || obj_name.compare(wl->this_wall->parent_object->last_name) == 0) {
          col = collide_wall_wrap(&start,&delta,wl->this_wall,&col_time,&col_hitpt,0);
          if ((col == COLLIDE_FRONT || col == COLLIDE_BACK) && col_time < t_sv_hit) {
            first_time = col_time;
            hit = col;
          } 
        }
      }
    }
    if (t_sv_hit != 1.0 && first_time < FOREVER)  {
      if (hit == COLLIDE_BACK) return 1;
      else return 0;
    }
  }
  return 0;
}

// Precession is handled using mcell_update_spins, but we want to keep track of
// magnetization to properly handle relaxation. This will update magnetization using
// the magnitude of Mxy but rotating to the match the computed precession angle. Should
// only be called after a diffusion gradient. Note that relaxation is independent of
// the precession, so we can update relaxation seperately at each time step or (I think) before
// any rf_pulses are applied.
int update_magnetization(mcellSim &sim, int nDirs)
{
  int i = 0, j, k, ns = 1;
  complex<double> I(0.0, 1.0);

  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  double Mxm, Mym, Mzm;
  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            if (sim.periodic || (fabs(0.01*mp->pos.x) < sim.xDim/2.0
                                 && fabs(0.01*mp->pos.y) < sim.yDim/2.0
                                 && fabs(0.01*mp->pos.z) < sim.zDim/2.0)) {
              k = mp->properties->species_id+1; //offset from total at 0
              for (j = 0; j < nDirs; j++) {
                Mxm = mp->Mx[j];
                Mym = mp->My[j];
                Mzm = mp->Mz[j];
                //mp->Mx[j] = (Mxm*cos(mp->phase_shift[j]) - Mym*sin(mp->phase_shift[j]));
                //mp->My[j] = (Mxm*sin(mp->phase_shift[j]) + Mym*cos(mp->phase_shift[j]));
                //mp->Mz[j] = Mzm;
                mp->Mx[j] = sqrt(Mxm*Mxm + Mym*Mym)*cos(mp->phase_shift[j]);
                mp->My[j] = sqrt(Mxm*Mxm + Mym*Mym)*sin(mp->phase_shift[j]);
                mp->Mz[j] = Mzm;
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

int update_relaxation(mcellSim &sim, int currentTime, int nDirs)
{
  int i = 0, j, ns = 1;
  complex<double> I(0.0, 1.0);

  /* loop over molecules and add signal from those in voxel */
  int activeMolecules = 0;
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  double Mxm, Mym, Mzm;
  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            if (sim.periodic || (fabs(0.01*mp->pos.x) < sim.xDim/2.0
                                 && fabs(0.01*mp->pos.y) < sim.yDim/2.0
                                 && fabs(0.01*mp->pos.z) < sim.zDim/2.0)) {
              activeMolecules++;
              for (j = 0; j < nDirs; j++) {
                mp->Mx[j] *= exp(-(double)currentTime/mp->properties->relax_t2);
                mp->My[j] *= exp(-(double)currentTime/mp->properties->relax_t2);
                mp->Mz[j] = mp->Mz[j]*exp(-(double)currentTime/mp->properties->relax_t1) + (1 - exp(-(double)currentTime/mp->properties->relax_t1));
              }
            }
          }
        }
      }
    }
  }
  return 0;
}

/* update position of all diffusing molecules */
void mcell_reset_origpos(mcellSim *s, int nDirs) {
                        
  int i = 0, j;
  struct ligand *curr_lig;
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;

  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            for (j = 0; j < nDirs; j++) {
              //mp->phase_shift[j] = 0.;
              mp->orig_pos.x = mp->pos.x;
              mp->orig_pos.y = mp->pos.y;
              mp->orig_pos.z = mp->pos.z;
            }
          }
        }
      }
    }
  }
}

int set_nc_current(mcellSim *s, input_data *d)
{
  mcellSim &sim = *(mcellSim *)s;
  int i, j, k, inside, nx = sim.ncNx, ny = sim.ncNy, nz = sim.ncNz;
  double x, y, z, dr, rx, ry, rz;
  ofstream ncjout; 

  if (sim.ncCurrent==0) return 0;
  sim.ncCurrentTimes.clear();
  sim.ncCurrentValues.clear();
  if (d->ncStrengthFile.size()>0) {
    ifstream ncFile(d->ncStrengthFile.c_str());
    double ts, nc;
    while (ncFile.good()) {
      ncFile >> ts >> nc;
      sim.ncCurrentTimes.push_back(ts);
      sim.ncCurrentValues.push_back(nc);
    } 
    ncFile.close();
  }
  if (sim.ncCurrentTimes.size()>0) sim.ncCurrent = 1;


  switch(sim.structurePacking) {
  case HEX2D:
    double dr,rx,ry,rz;
    sim.vp = new CVectorPotential(nx,ny,1,sim.xDim,sim.yDim,sim.zDim);
    z = 0;
    for (j=0; j<ny; j++) {
      y = -sim.yDim/2.0 + j*sim.yDim/(double)ny;
      for (i=0; i<nx; i++) {
        x = -sim.xDim/2.0 + i*sim.xDim/(double)nx;
        rx = ry = 0;
        k = 0;
        dr = (x-rx)*(x-rx)+(y-ry)*(y-ry);
        if (sqrt(dr)<sim.cylinderRadius) k=1;
        if (k==0) {
          rx = -sim.xDim/2; ry = -sim.yDim/2;
          dr = (x-rx)*(x-rx)+(y-ry)*(y-ry);
          if (sqrt(dr)<sim.cylinderRadius) k=1;
        }
        if (k==0) {
          rx = sim.xDim/2; ry = -sim.yDim/2;
          dr = (x-rx)*(x-rx)+(y-ry)*(y-ry);
          if (sqrt(dr)<sim.cylinderRadius) k=1;
        }
        if (k==0) {
          rx = -sim.xDim/2; ry = sim.yDim/2;
          dr = (x-rx)*(x-rx)+(y-ry)*(y-ry);
          if (sqrt(dr)<sim.cylinderRadius) k=1;
        }
        if (k==0) {
          rx = sim.xDim/2; ry = sim.yDim/2;
          dr = (x-rx)*(x-rx)+(y-ry)*(y-ry);
          if (sqrt(dr)<sim.cylinderRadius) k=1;
        }
        if (k==0) sim.vp->getJ()[i+nx*j] = 1.0;
        else sim.vp->getJ()[i+nx*j] = 0;
      }
    } 
    sim.vp->solveA();
    sim.vp->printBm(0, d->ncOutputFile.c_str());
    break;
  case ANYREGION:
    sim.vp = new CVectorPotential(nx,ny,nz,sim.xDim,sim.yDim,sim.zDim, d->BCx0, d->BCx1, d->BCy0, d->BCy1, d->BCz0, d->BCz1);
    if (d->ncJinFile.size()>0) {
      cout << "Loading ncCurrent field" << endl;
      ifstream ncjin(d->ncJinFile.c_str());
      int position = 0;
      double curval = 0.;
      while (ncjin.good()) {
        ncjin >> curval;
        sim.vp->getJ()[position] = curval;
        position++;
      }
      ncjin.close();
    } else {
      if (d->ncJoutFile.size()>0) ncjout.open(d->ncJoutFile.c_str());
      for (k=0; k<nz; k++) {
        z = -sim.zDim/2.0 + k*sim.zDim/(double)nz;
        cout << "k coord: " << k << endl;
        for (j=0; j<ny; j++) {
          y = -sim.yDim/2.0 + j*sim.yDim/(double)ny;
          for (i=0; i<nx; i++) {
            x = -sim.xDim/2.0 + i*sim.xDim/(double)nx;
            inside = 0;
            if (test_point_inside(x,y,z, d->ncObjectName, s)) inside = 1;
            if (inside==0) {
              sim.vp->getJ()[k+nz*(j+ny*i)] = 1.0;
              if (d->ncJoutFile.size()>0) ncjout << 1. << endl;
            } else {
              sim.vp->getJ()[k+nz*(j+ny*i)] = 0;
              if (d->ncJoutFile.size()>0) ncjout << 0. << endl;
            }
          }
        }
      }
      if (d->ncJoutFile.size()>0) ncjout.close();
    }
    sim.vp->solveA();
    sim.vp->printBm(0, d->ncOutputFile.c_str());
    break;

  default:
    sim.ncCurrent=0;
    return 0;
  }
  return 1;
}

extern "C" {
  int is_in_region(struct volume_molecule *mp)
  {
    mcellSim &sim = *(mcellSim *)world->msim;
    double dr,rx,ry,rz;
    double cx[14],cy[14],cz[14];
    int i, nc;
    int k;
    double xs = sim.xScale, ys = sim.yScale, zs = sim.zScale;
    double sx=sim.xDim/2.0*100, sy=sim.yDim/2.0*100, sz=sim.zDim/2.0*100;
    if ((strcmp(mp->properties->sym->name,"water_internal")!=0)&&
        (strcmp(mp->properties->sym->name,"water_external")!=0))
      return 1;
    switch(sim.structurePacking) {
    case NONE: return 1;
    case HEX2D:
      rx = ry = 0;
      k = 0; // 0 - external; 1 - internal
      dr = (mp->pos.x-rx)*(mp->pos.x-rx)+(mp->pos.y-ry)*(mp->pos.y-ry);
      if (sqrt(dr)<sim.cylinderRadius*100) k=1;
#if 1
      if (k==0) {
        rx = -sim.xDim/2*100; ry = -sim.yDim/2*100;
        dr = (mp->pos.x-rx)*(mp->pos.x-rx)+(mp->pos.y-ry)*(mp->pos.y-ry);
        if (sqrt(dr)<sim.cylinderRadius*100) k=1;
      }
      if (k==0) {
        rx = sim.xDim/2*100; ry = -sim.yDim/2*100;
        dr = (mp->pos.x-rx)*(mp->pos.x-rx)+(mp->pos.y-ry)*(mp->pos.y-ry);
        if (sqrt(dr)<sim.cylinderRadius*100) k=1;
      }
      if (k==0) {
        rx = -sim.xDim/2*100; ry = sim.yDim/2*100;
        dr = (mp->pos.x-rx)*(mp->pos.x-rx)+(mp->pos.y-ry)*(mp->pos.y-ry);
        if (sqrt(dr)<sim.cylinderRadius*100) k=1;
      }
      if (k==0) {
        rx = sim.xDim/2*100; ry = sim.yDim/2*100;
        dr = (mp->pos.x-rx)*(mp->pos.x-rx)+(mp->pos.y-ry)*(mp->pos.y-ry);
        if (sqrt(dr)<sim.cylinderRadius*100) k=1;
      }
#endif
      if ((strcmp(mp->properties->sym->name,"water_internal")==0)&&(k==1)) 
        return 1;
      if ((strcmp(mp->properties->sym->name,"water_external")==0)&&(k==0)) 
        return 1;
      return 0;
    case CCP3D:
      k = 0;  cx[k] = 0; cy[k] = 0; cz[k] = -sz;
      k = 1;  cx[k] = 0; cy[k] = 0; cz[k] =  sz;
      k = 2;  cx[k] = 0; cy[k] = -sy; cz[k] = 0;
      k = 3;  cx[k] = 0; cy[k] =  sy; cz[k] = 0;
      k = 4;  cx[k] = -sx; cy[k] = 0; cz[k] = 0;
      k = 5;  cx[k] =  sx; cy[k] = 0; cz[k] = 0;
      k = 6;  cx[k] = -sx; cy[k] = -sy; cz[k] = -sz;
      k = 7;  cx[k] = -sx; cy[k] = -sy; cz[k] =  sz;
      k = 8;  cx[k] = -sx; cy[k] =  sy; cz[k] = -sz;
      k = 9;  cx[k] = -sx; cy[k] =  sy; cz[k] =  sz;
      k = 10; cx[k] =  sx; cy[k] = -sy; cz[k] = -sz;
      k = 11; cx[k] =  sx; cy[k] = -sy; cz[k] =  sz;
      k = 12; cx[k] =  sx; cy[k] =  sy; cz[k] = -sz;
      k = 13; cx[k] =  sx; cy[k] =  sy; cz[k] =  sz;
      nc = k+1;
      //k = 0;  cx[k] = 0; cy[k] = 0; cz[k] = 0;
      //nc = 1;
      for (i = 0; i<nc; i++) {
        k = 0; // 0 - external; 1 - internal
        rx = cx[i]; ry = cy[i]; rz = cz[i];
        dr = (mp->pos.x-rx)*(mp->pos.x-rx)/xs/xs+
             (mp->pos.y-ry)*(mp->pos.y-ry)/ys/ys+
             (mp->pos.z-rz)*(mp->pos.z-rz)/zs/zs;
        if (sqrt(dr)<sim.sphereRadius*100) k=1;
        if (k==1) break;
      }
      if ((strcmp(mp->properties->sym->name,"water_internal")==0)&&(k==1)) 
        return 1;
      if ((strcmp(mp->properties->sym->name,"water_external")==0)&&(k==0)) 
        return 1;
      return 0;
    case HCP3D:
      return 1;
    }
  }
}

void count_regions(int iteration)
{
  int i = 0, j, k, k0, ns = 1, n_species, nc;
  int activeMolecules = 0;
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  std::vector<int> ind;
  std::vector<int> cnt;
  std::vector<int> cnt1;
  ind.assign(world->n_species,0);
  mcellSim &sim = *(mcellSim *)world->msim;

#if 1
  for (j=0; j<world->n_species; j++) 
    if (world->species_list[j]->population) {
      ind[j]=ns;
#ifdef DEBUG
      cout << " " << ns << " " << j << " " 
           << world->species_list[j]->population << endl;
#endif
      ns++;
    }
#endif
  if (ns<3) ns = 3;
  n_species = ns;
  cnt.assign(ns,0);
  cnt1.assign(ns,0);
  nc = 0;
  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            if (sim.periodic || (fabs(0.01*mp->pos.x) < sim.xDim/2.0
                                 && fabs(0.01*mp->pos.y) < sim.yDim/2.0
                                 && fabs(0.01*mp->pos.z) < sim.zDim/2.0)) {
              activeMolecules++;
              k = ind[mp->properties->species_id];
              cnt[k]++;
              k0 = k;
              if (strcmp(mp->properties->sym->name,"water_internal")==0) {
                if (is_in_region(mp)) k = 1; else k = 2;
              } else {
                if (is_in_region(mp)) k = 2; else k = 1;
              }
#if 1
              if (k==0) k=2;
              cnt1[k]++;
              if (k!=k0) {
                nc++;
#ifdef DEBUG
                std::cout << "@ " << k << "<>" << k0 << " "
                          << mp->pos.x/100 << " " << mp->pos.y/100 << " "
                          << mp->pos.z/100 << " " 
                          << mp->orig_pos.x/100 << " " 
                          << mp->orig_pos.y/100 << " "
                          << mp->orig_pos.z/100 << " " 
                          << mp->properties->sym->name << " "
                          << mp << endl;
#endif
              }
#endif
            }
          }
        }
      }
    }
  }
#ifdef DEBUG
  for (j=1; j<ns; j++) cnt[0] += cnt[j];
  for (j=1; j<ns; j++) cnt1[0] += cnt1[j];
  cout << "active molecules: " << activeMolecules << "/" << nc;
  for (j=0; j<ns; j++) cout << " " << cnt[j] << " (" << cnt1[j] << ")";
  cout << endl;
#endif
}

/* compute signal from particle phase shifts -- returns 0 if successful */
int Compute_Signal(complex<double> *&signal, // complex signal 
                   double *&VSignal,   // vertex signal
                   mcellSim &sim,
                   int num_gradient_directions,
                   double total_weight,
                   double g_snr, int rank, int n_species) {

  int i = 0, j, k, ns = 1;
  std::vector<int> cnt;
  complex<double> I(0.0, 1.0);
 
  ns = n_species;
  cnt.assign(ns,0);

  /* start with zero signal */
  signal = new complex<double>[num_gradient_directions*ns];
  VSignal = new double[num_gradient_directions*ns];
  for (j = 0; j < num_gradient_directions*ns; j++) {
    signal[j] = 0.0;
  }

  /* loop over molecules and add signal from those in voxel */
  int activeMolecules = 0;
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct volume_molecule *mp;
  for (slp = world->storage_head ; slp != NULL ; slp = slp->next) {
    for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale) {
      for (i=-1;i<shp->buf_len;i++) {
        for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next) {
          amp = (struct abstract_molecule*)aep;
          if (amp->properties == NULL) continue;
          if ((amp->properties->flags & NOT_FREE)==0) {
            mp = (struct volume_molecule*)amp;
            if (sim.periodic || (fabs(0.01*mp->pos.x) < sim.xDim/2.0
                                 && fabs(0.01*mp->pos.y) < sim.yDim/2.0
                                 && fabs(0.01*mp->pos.z) < sim.zDim/2.0)) {
              activeMolecules++;
              k = mp->properties->species_id+1; //offset from total held at index 0
              cnt[k]++;

              for (j = 0; j < num_gradient_directions; j++) {
                // Must be a way to do this cleaner without tempmag, whats syntax?
                complex<double> tempmag(mp->Mx[j], mp->My[j]);
                signal[j] += tempmag;
                signal[j+k*num_gradient_directions] += tempmag;
              }
            }
          }
        }
      }
    }
  }
  
  for (j=1; j<ns; j++) cnt[0] += cnt[j];
  /* scaling for something <= 1 */
  double renorm;
  for (k = 0; k < ns; k++) {
    if (cnt[k]<1) continue;
    renorm = 1.0 / (double)cnt[k] / total_weight;
    for (j = 0; j < num_gradient_directions; j++) {
      signal[j+k*num_gradient_directions] *= renorm;
    }
  }
  
  for (int i=0; i<num_gradient_directions*ns; i++) { // for each spin ...
    VSignal[i] = abs(signal[i]);
  }
  return 0;
}

float CurrentGradientStrength(int time, sequenceTimes sTimes, mcellSim sim) {
  float G = sim.gradient_strength;
  float cg = 0.0;
  if (time >= sTimes.rampUp1 && time < sTimes.gOn1) {
    /* first up ramp */
    cg = G *
      ((double)time / (double)sim.gradient_ramp_time);
  }
  if (time >= sTimes.rampUp2 && time < sTimes.gOn2) {
    /* second up ramp */
    cg = G *
      ((double)(time - sTimes.rampUp2)
       / (double)sim.gradient_ramp_time);
  }
  if (time >= sTimes.rampDown1 && time < sTimes.gOff) {
    /* first down ramp */
    cg = G*
      (1.0 - ((double)(time - sTimes.rampDown1)
              / (double)sim.gradient_ramp_time));
  }
  if (time >= sTimes.rampDown2 && time < sTimes.end) {
    /* second down ramp */
    cg = G*
      (1.0 - ((double)(time - sTimes.rampDown2)
              / (double)sim.gradient_ramp_time));
  }
  if ((time >= sTimes.gOn1 && time < sTimes.rampDown1)
      || (time >= sTimes.gOn2 && time < sTimes.rampDown2)) {
    /* full strength gradient */
    cg = G;
  }
  return cg;
}

int LoadTensorFile(input_data *d, float *gx, float *gy, float *gz, int nDirs) {

    ifstream qf;
    istringstream istr;
    string s0;
    int i, nvec;
    if (nDirs<=0) return 0;
    qf.open(d->tensor_file.c_str());
    if (!qf.is_open()) {
      cerr << " Error openning tensor file '" << d->tensor_file << endl;
      exit(-1);
    }

    while (qf.good()) {
      s0.clear();
      getline(qf, s0);
      if (!qf.good()) break;
      if (s0[0]=='#') continue;
      istr.str(s0);
      istr >> nvec;
      istr.clear();
      if  (nDirs != nvec) {
        for (i=0; i<nvec; i++) getline(qf, s0);
      } else {
        for (i=0; i<nvec; i++) {
          s0.clear();
          getline(qf, s0);
          istr.str(s0);
          istr >> gx[i] >> gy[i] >> gz[i];
          istr.clear();
        }
        qf.close();
        return nvec;
      }
    }
    qf.close();
    return -1;
}

int SaveSignalFile(mcellSim &sim, input_data *d, int nDirs, int currentTime) {
    int rank = 0;
    complex<double> *complexSignal;
    double *vertex_signal;
    int n_species = world->n_species + 1; //0 is for total signal
    if(Compute_Signal(complexSignal, vertex_signal, sim, nDirs, 1.0, d->snr, rank, n_species)) {
      cerr << " Error computing signal " << endl;
      exit(-1);
    }
    
    struct species *sp;  
    if (rank == 0) {
      if (d->complex_signal_filename.length()) {
        /* write out complex signal */
        ofstream complexSignalFile(d->complex_signal_filename.c_str(), ios::out | ios::app);
        complexSignalFile.precision(12);
        complexSignalFile << "# Signal at " << currentTime << " us \n";
        for (int i = 0; i < nDirs; i++) {
          if (i == 0) { //Write labels and b0
            complexSignalFile << "# ";
            for (int k = 0; k < n_species; k++) {
              if (k == 0) {
                complexSignalFile << "total ";
              } else {
                sp = world->species_list[k-1]; 
                if (sp != world->all_mols && sp != world->all_volume_mols && sp != world->all_surface_mols) {
                  complexSignalFile <<  world->species_list[k-1]->sym->name << " ";
                }
              }
            }
            complexSignalFile << endl;
          }
          for (int k = 0; k < n_species; k++) {
            if (k == 0) { 
              complexSignalFile << real(complexSignal[i+k*nDirs]) << " " 
                                << imag(complexSignal[i+k*nDirs]) << " ";
            } else {
              sp = world->species_list[k-1];
              if (sp->population && sp != world->all_mols && sp != world->all_volume_mols && sp != world->all_surface_mols) {
                complexSignalFile << real(complexSignal[i+k*nDirs]) << " " 
                                << imag(complexSignal[i+k*nDirs]) << " ";
              }
            }
          }
          complexSignalFile << endl;
        }
        complexSignalFile.close();
      }
      if (d->signal_filename.length()) {
        /* write out the vertex signal */
        ofstream vertex_signal_file(d->signal_filename.c_str(), ios::out | ios::app);
        vertex_signal_file.precision(12);
        vertex_signal_file << "# Signal at " << currentTime << " us\n";
        for (int i = 0; i < nDirs; i++) {
          if (i == 0) {
            vertex_signal_file << "# ";
            for (int k = 0; k < n_species; k++) {
              if (k == 0) {
                vertex_signal_file << "total" << " ";
              } else {
                sp = world->species_list[k-1];
                if (sp->population && sp != world->all_mols && sp != world->all_volume_mols && sp != world->all_surface_mols) {
                  vertex_signal_file <<  world->species_list[k-1]->sym->name << " ";
                }
              }
            }
            vertex_signal_file << endl;
          }
          for (int k = 0; k < n_species; k++) {
            if (k == 0 ) { 
              vertex_signal_file << vertex_signal[i+k*nDirs] << " ";
            } else {
              sp = world->species_list[k-1];
              if (sp->population && sp != world->all_mols && sp != world->all_volume_mols && sp != world->all_surface_mols) {
                vertex_signal_file << vertex_signal[i+k*nDirs] << " ";
              }
            }
          }
          vertex_signal_file << endl;
        }
        vertex_signal_file.close();
      }
    }
  return 0;
}

void set_sequenceTimes(mcellSim &sim, sequenceTimes *sTimes) {
  int gradient_time, echo_remainder;

  switch (sim.pulse_sequence) {
  case PULSE_SPIN_ECHO:
  case PULSE_GRAD_ECHO:
  case PULSE_STEAM_DTI:
    //sTimes = new sequenceTimes[1];
    sTimes[0].rampUp1 = 0;
    sTimes[0].gOn1 = sTimes[0].rampUp1 + sim.gradient_ramp_time;
    sTimes[0].rampDown1 = sTimes[0].gOn1 + sim.gradient_on_time;
    sTimes[0].gOff = sTimes[0].rampDown1 + sim.gradient_ramp_time;
    sTimes[0].rampUp2 = sTimes[0].gOff + sim.gradient_off_time;
    sTimes[0].gOn2 = sTimes[0].rampUp2 + sim.gradient_ramp_time;
    sTimes[0].rampDown2 = sTimes[0].gOn2 + sim.gradient_on_time;
    sTimes[0].end = sTimes[0].rampDown2 + sim.gradient_ramp_time;
    sim.endTime = sTimes[0].end;
    if (sim.repeat_time < sim.endTime) sim.repeat_time = sim.endTime;
    if (sim.pulse_sequence == PULSE_STEAM_DTI) {
      if (sim.mixing_time <= 1) {
        cout << "WARN: STEAM DTI mixing time is zero .... changed to SPIN ECHO"
             << endl;
        sim.pulse_sequence = PULSE_SPIN_ECHO;
        sim.mixing_time = 0;
      }
      if (sim.mixing_time > sim.gradient_off_time)
        sim.mixing_time = sim.gradient_off_time;
    }
    break;
  case PULSE_MULTI_SPIN_ECHO:
    gradient_time = 2*sim.gradient_on_time + sim.gradient_off_time + 4*sim.gradient_ramp_time;
    if (sim.echo_time > gradient_time) {
      echo_remainder = (sim.echo_time - gradient_time)/2;
    } else {
      echo_remainder = 0; 
      sim.echo_time = gradient_time;
    }
    sTimes[0].rampUp1 = echo_remainder;
    sTimes[0].gOn1 = sTimes[0].rampUp1 + sim.gradient_ramp_time;
    sTimes[0].rampDown1 = sTimes[0].gOn1 + sim.gradient_on_time;
    sTimes[0].gOff = sTimes[0].rampDown1 + sim.gradient_ramp_time;
    sTimes[0].rampUp2 = sTimes[0].gOff + sim.gradient_off_time;
    sTimes[0].gOn2 = sTimes[0].rampUp2 + sim.gradient_ramp_time;
    sTimes[0].rampDown2 = sTimes[0].gOn2 + sim.gradient_on_time;
    sTimes[0].end = sTimes[0].rampDown2 + sim.gradient_ramp_time;
    sim.endTime = sTimes[0].end + echo_remainder;
    break;
  case PULSE_DPFG:
    //sTimes = new sequenceTimes[2];
    sTimes[0].rampUp1 = 0;
    sTimes[0].gOn1 = sTimes[0].rampUp1 + sim.gradient_ramp_time;
    sTimes[0].rampDown1 = sTimes[0].gOn1 + sim.gradient_on_time;
    sTimes[0].gOff = sTimes[0].rampDown1 + sim.gradient_ramp_time;
    sTimes[0].rampUp2 = sTimes[0].gOff + sim.gradient_off_time;
    sTimes[0].gOn2 = sTimes[0].rampUp2 + sim.gradient_ramp_time;
    sTimes[0].rampDown2 = sTimes[0].gOn2 + sim.gradient_on_time;
    sTimes[0].end = sTimes[0].rampDown2 + sim.gradient_ramp_time;
    sTimes[1].rampUp1 = sTimes[0].rampUp2 + sim.mixing_time;
    sTimes[1].gOn1 = sTimes[1].rampUp1 + sim.gradient_ramp_time;
    sTimes[1].rampDown1 = sTimes[1].gOn1 + sim.gradient_on_time;
    sTimes[1].gOff = sTimes[1].rampDown1 + sim.gradient_ramp_time;
    sTimes[1].rampUp2 = sTimes[1].gOff + sim.gradient_off_time;
    sTimes[1].gOn2 = sTimes[1].rampUp2 + sim.gradient_ramp_time;
    sTimes[1].rampDown2 = sTimes[1].gOn2 + sim.gradient_on_time;
    sTimes[1].end = sTimes[1].rampDown2 + sim.gradient_ramp_time;
    sim.endTime = sTimes[1].end;
    if (sim.repeat_time < sim.endTime) sim.repeat_time = sim.endTime;
    break;
  case PULSE_DWSSFP:
    //sTimes = new sequenceTimes[1];
    sTimes[0].rampUp1 = 0;
    sTimes[0].gOn1 = sTimes[0].rampUp1 + sim.gradient_ramp_time;
    sTimes[0].rampDown1 = sTimes[0].gOn1 + sim.gradient_on_time;
    sTimes[0].gOff = sTimes[0].rampDown1 + sim.gradient_ramp_time;
    sTimes[0].rampUp2 = sTimes[0].gOff;
    sTimes[0].gOn2 = sTimes[0].gOff;
    sTimes[0].rampDown2 = sTimes[0].gOff;
    sTimes[0].end = sTimes[0].gOff;
    sim.endTime = sTimes[0].end;
    if (sim.repeat_time < sim.endTime) sim.repeat_time = sim.endTime;
    break;
  }
}

void simulation(input_data *d) {
  int rank = 0;
  int nProcessors = 1;
  int iteration = 0;
  int NSTEPS;//, step_size;
  double step_size;
  int nDirs;
  int currentTime, currentSeqTime; //is time in the current sequence if there are multiple repeats
  mcellSim sim;
  double rand_prec;

  complex<double> *complexSignal;
  double *vertex_signal;
  double *gvec_xaxis;
  
  ostringstream os;
  double maxDx, maxDy, maxDz;
  int debug = 0;
  float *gx, *gy, *gz;
  int nparticles;
  double (*gradient_directions)[3];
  double (*dpfg_gradient_directions)[3]; // Holds actual gradients for DPFG, where gradient_directions holds
                                         // the vector sum of these gradients and a gradient along the x axis

  time_t t_start = time(NULL), t_end, t_now;
  sim.debugMCell = d->debugMCell;
  sim.gradient_off_time = d->Delta;
  sim.gradient_on_time = d->delta;
  sim.current_gradient_strength = d->gradient_strength;
  sim.gradient_strength = d->gradient_strength;
  sim.bath_concentration = d->bath_concentration;
  sim.sheath_concentration = d->sheath_concentration;
  sim.core_concentration = d->core_concentration;
  sim.D_bath = d->D_bath;
  sim.D_sheath = d->D_sheath;
  sim.D_core = d->D_core;
  sim.gradient_ramp_time = d->ramp;
  sim.pulse_sequence = d->pulse;
  sim.mixing_time = d->mixing_time;
  sim.dpfgAngle = d->dpfgAngle;
  sim.EulerAngle1 = d->EulerAngle1;
  sim.EulerAngle2 = d->EulerAngle2;
  sim.EulerAngle3 = d->EulerAngle3;
  sim.ncEulerAngle1 = d->ncEulerAngle1;
  sim.ncEulerAngle2 = d->ncEulerAngle2;
  sim.ncEulerAngle3 = d->ncEulerAngle3;
  sim.periodic = d->periodic;
  sim.parallelPlates = d->parallelPlates;
  sim.voxel_size = d->voxel_size;
  sim.volume_size = d->volume_size;
  sim.xDim = d->xDim;
  sim.yDim = d->yDim;
  sim.zDim = d->zDim;
  sim.xScale = d->xScale;
  sim.yScale = d->yScale;
  sim.zScale = d->zScale;
  sim.volumeHeight = d->volumeHeight;
  sim.coreRadius = d->coreRadius;
  sim.sheathRadius = d->sheathRadius;
  sim.sphereRadius = d->sphereRadius;
  sim.cylinderRadius = d->cylinderRadius;
  sim.cylinderSpacing = d->cylinderSpacing;
  sim.cylinderHeight = d->cylinderHeight;
  sim.threeRegion = d->threeRegion;
  sim.permeable = d->permeable;
  sim.P_scale = d->P_scale;;
  sim.ncCurrent = d->ncCurrent;;
  sim.cylinderArray = d->cylinderArray;
  sim.structureArray = d->structureArray;
  sim.structurePacking = d->structurePacking;
  NSTEPS = d->nsteps;
  sim.timestep = d->step_size;
  sim.orig_timestep = d->step_size;
  sim.num_particles = d->number_of_particles;
  sim.gradientDirections = d->gradientDirections;
  sim.ncNx = d->nc_Nx;
  sim.ncNy = d->nc_Ny;
  sim.ncNz = d->nc_Nz;
  sim.use_relaxation = d->use_relax;
  sim.repeat_number = d->repeat_number;
  sim.repeat_time = d->repeat_time;
  sim.flip_angle = d->flip_angle;
  sim.echo_time = d->echo_time;
  sim.echo_spacing = d->echo_spacing;
  sim.echo_number = d->echo_number;

  if (d->mcell_file == "" && d->structureFile == "") {
    cout << "No MDL file or structure file specified: calculating diffusion within\
bounding box.." << endl;
  }
  if (rank == 0) {
    printf("\nGenerating gradients: ... ");
  }
  if (sim.pulse_sequence == PULSE_DPFG) {
    /* BMR - gradient_directions now holds a sum of the gradient along the x axis and 
     * the other specified gradients, which will now be kept in a separate array. This
     * allows multiple gradients to be run in a single simulation
     *
    /* PULSE_DPFG is a special case: there is only one gradient direction
       at a time---so nDirs = 1---but it may be a composite when gradients
       overlap. */
    sim.dpfgAngleList.clear();
    if (d->dpfgAngleFile.size()>0) {
      ifstream dpfgFile(d->dpfgAngleFile.c_str());
      double angle;
      while (dpfgFile.good()) {
        dpfgFile >> angle;
        sim.dpfgAngleList.push_back(angle);
      } 
      dpfgFile.close();
      if (nDirs != sim.dpfgAngleList.size()) {
        nDirs = sim.dpfgAngleList.size(); 
      }
      gx = new float[nDirs];
      gy = new float[nDirs];
      gz = new float[nDirs];
      for (int dir = 0; dir < nDirs; dir++) {
      /* These hold the directions for each gradient echo pair */
        gx[dir] = cos(M_PI*sim.dpfgAngleList[dir]/180.0);
        gy[dir] = sin(M_PI*sim.dpfgAngleList[dir]/180.0);
        gz[dir] = 0.0;
      }
    } else if (!d->tensor_file.empty() && d->tensorFileDPFG) { //Specify both vectors in tensor file 
      cout << "Using DPFG tensor file";
      ifstream dpfgFile(d->tensor_file.c_str());
      double gx1,gy1,gz1,gx2,gy2,gz2;
      sim.gx1.clear(); sim.gy1.clear(); sim.gz1.clear();
      sim.gx2.clear(); sim.gy2.clear(); sim.gz2.clear();
      //while (dpfgFile.good()) {
      while(dpfgFile >> gx1 >> gy1 >> gz1 >> gx2 >> gy2 >> gz2) {
        sim.gx1.push_back(gx1);
        sim.gy1.push_back(gy1);
        sim.gz1.push_back(gz1);
        sim.gx2.push_back(gx2);
        sim.gy2.push_back(gy2);
        sim.gz2.push_back(gz2);
      } 
      dpfgFile.close();
      if (nDirs != sim.gx1.size()) {
        nDirs = sim.gx1.size(); 
      }
      // To try to protect from bad behavior, make second vector fill gx,gy,gz, but really want
      // to later use sim.g** for later computation
      gx = new float[nDirs];
      gy = new float[nDirs];
      gz = new float[nDirs];
      for (int dir = 0; dir < nDirs; dir++) {
      /* These hold the directions for each gradient echo pair */
        gx[dir] = sim.gx2[dir];
        gy[dir] = sim.gy2[dir];
        gz[dir] = sim.gz2[dir];
      }
    } else if (!d->tensor_file.empty()) {
      nDirs = sim.gradientDirections;
      gx = new float[nDirs];
      gy = new float[nDirs];
      gz = new float[nDirs];
      if (LoadTensorFile(d, gx, gy, gz, nDirs) != nDirs) {
        cerr << " Error loading tensor file '" << d->tensor_file << endl;
        exit(-1);
      }
    } else {
      sim.gradientDirections = 1;
      nDirs = 1;
      gx = new float[nDirs];
      gy = new float[nDirs];
      gz = new float[nDirs];
      gx[0] = cos(M_PI*sim.dpfgAngle/180.0);
      gy[0] = sin(M_PI*sim.dpfgAngle/180.0);
      gz[0] = 0.0;
    } 
  }
  else if (!d->tensor_file.empty()) {
    nDirs = sim.gradientDirections;
    gx = new float[nDirs];
    gy = new float[nDirs];
    gz = new float[nDirs];
    if (LoadTensorFile(d, gx, gy, gz, nDirs) != nDirs) {
      cerr << " Error loading tensor file '" << d->tensor_file << endl;
      exit(-1);
    }
  }
  else if (d->ntess > 0) {
    nDirs = 5*(1 << (2*d->ntess - 1)) + 2;
    int ntris  = 20*(1 << (2*d->ntess));
    int *tris = new int[ntris];
    gx = new float[nDirs];
    gy = new float[nDirs];
    gz = new float[nDirs];
    assert(sphtes(gx, gy, gz, tris, POLY_TYPE, d->ntess)
           == nDirs);
    delete [] tris;
  }
  else {
    /* around a semicircle in the x-z plane */
    nDirs = sim.gradientDirections;
    gx = new float[nDirs];
    gy = new float[nDirs];
    gz = new float[nDirs];
    int dir;
    for (dir = 0; dir < nDirs; dir++) {
      gx[dir] = cos((M_PI*dir)/sim.gradientDirections);
      gy[dir] = 0.0;
      gz[dir] = sin((M_PI*dir)/sim.gradientDirections);
      cout << "Gradient " << dir + 1 
           << ": [" << gx[dir] << ", " 
           << gy[dir] << ", " 
           << gz[dir] << "]" << endl;
    }
    // switch(sim.gradientDirections) 
    // case 1:
    //   nDirs = 1;
    //   gx = new float[nDirs];
    //   gy = new float[nDirs];
    //   gz = new float[nDirs];
    //   gx[0] = 1.0;
    //   gy[0] = 0.0;
    //   gz[0] = 0.0;
    //   break;
    // case 2:
    //   nDirs = 2;
    //   gx = new float[nDirs];
    //   gy = new float[nDirs];
    //   gz = new float[nDirs];
    //   gx[0] = 1.0;
    //   gy[0] = 0.0;
    //   gz[0] = 0.0;
    //   gx[1] = 0.0;
    //   gy[1] = 0.0;
    //   gz[1] = 1.0;
    //   break;
    // default:
    //   cerr << "Except when generating gradient directions using a tessellated icosahedron,\nthe number of gradient directions must be 1 or 2. Please choose either\nntess = <tessellation level> or \ngradient directions = <1 or 2>\n in your input file." << endl;
    //   exit(-1);
  }
#if 0
  if (rank == 0) {
    printf(" using %i gradient directions.\n", nDirs);
    if (!nDirs) {
      cerr << "Zero gradient directions.  Something went wrong." << endl;
      exit(-2);
    }
  }
#endif
  
  int gradientDirSize = nDirs;
  // Apply rotations with Euler angles
  float tx, ty, tz, a1, a2, a3;
  a1 = M_PI*sim.ncEulerAngle1/180.0;
  a2 = M_PI*sim.ncEulerAngle2/180.0;
  a3 = M_PI*sim.ncEulerAngle3/180.0;
  sim.ncTr[0][0] = cos(a1) * cos(a3) - sin(a1) * cos(a2) * sin(a3);
  sim.ncTr[0][1] = sin(a1) * cos(a3) + cos(a1) * cos(a2) * sin(a3);
  sim.ncTr[0][2] = sin(a2) * sin(a3);
  sim.ncTr[1][0] = -cos(a1) * sin(a3) - sin(a1) * cos(a2) * cos(a3);
  sim.ncTr[1][1] = -sin(a1) * sin(a3) + cos(a1) * cos(a2) * cos(a3);
  sim.ncTr[1][2] = sin(a2) * cos(a3);
  sim.ncTr[2][0] = sin(a1) * sin(a2);
  sim.ncTr[2][1] = -cos(a1) * sin(a2);
  sim.ncTr[2][2] = cos(a2);
  a1 = M_PI*sim.EulerAngle1/180.0;
  a2 = M_PI*sim.EulerAngle2/180.0;
  a3 = M_PI*sim.EulerAngle3/180.0;
  sim.Tr[0][0] = cos(a1) * cos(a3) - sin(a1) * cos(a2) * sin(a3);
  sim.Tr[0][1] = sin(a1) * cos(a3) + cos(a1) * cos(a2) * sin(a3);
  sim.Tr[0][2] = sin(a2) * sin(a3);
  sim.Tr[1][0] = -cos(a1) * sin(a3) - sin(a1) * cos(a2) * cos(a3);
  sim.Tr[1][1] = -sin(a1) * sin(a3) + cos(a1) * cos(a2) * cos(a3);
  sim.Tr[1][2] = sin(a2) * cos(a3);
  sim.Tr[2][0] = sin(a1) * sin(a2);
  sim.Tr[2][1] = -cos(a1) * sin(a2);
  sim.Tr[2][2] = cos(a2);
  for (int i = 0; i < gradientDirSize; i++) {
    tx = gx[i];
    ty = gy[i];
    tz = gz[i];
    gx[i] = sim.Tr[0][0]*tx + sim.Tr[0][1]*ty + sim.Tr[0][2]*tz;
    gy[i] = sim.Tr[1][0]*tx + sim.Tr[1][1]*ty + sim.Tr[1][2]*tz;
    gz[i] = sim.Tr[2][0]*tx + sim.Tr[2][1]*ty + sim.Tr[2][2]*tz;
  }
  gradient_directions = new double[gradientDirSize+1][3];
  
  for (int i = 0; i < gradientDirSize; i++) {
    if (rank == 0) {
      if(debug) printf("g[%i]=(%f,%f,%f)\n",i,gx[i],gy[i],gz[i]) ;
    }
    gradient_directions[i+1][0] = gx[i];
    gradient_directions[i+1][1] = gy[i];
    gradient_directions[i+1][2] = gz[i];
  }
  gradient_directions[0][0] = gradient_directions[0][1] =
    gradient_directions[0][2] = 0;
  if (sim.pulse_sequence == PULSE_DPFG) {
    dpfg_gradient_directions = new double[gradientDirSize+1][3];
    for (int i = 0; i < gradientDirSize; i++) {
      dpfg_gradient_directions[i][0] = gx[i];
      dpfg_gradient_directions[i][1] = gy[i];
      dpfg_gradient_directions[i][2] = gz[i];
    } 
    dpfg_gradient_directions[0][0] = dpfg_gradient_directions[0][1] =
      dpfg_gradient_directions[0][2] = 0;
  }
  gradientDirSize++;
  nDirs = sim.gradientDirections = gradientDirSize;


  delete [] gx;
  delete [] gy;
  delete [] gz;

  //Num particles not used anymore, should be specifying molecules through MDL to have different populations
  /*if (rank == 0) {
    cout << "Using " << sim.num_particles << " diffusing particles." << endl << endl;
  }*/

  /*printf("ncCurrent rotation matrix: [%0.3f,%0.3f,%0.3f;\n",sim.ncTr[0][0],sim.ncTr[0][1],sim.ncTr[0][2]);
  printf("                            %0.3f,%0.3f,%0.3f;\n",sim.ncTr[1][0],sim.ncTr[1][1],sim.ncTr[1][2]);
  printf("                            %0.3f,%0.3f,%0.3f]\n",sim.ncTr[2][0],sim.ncTr[2][1],sim.ncTr[2][2]);*/

  currentTime = 0;
  B0.t = currentTime;
  sequenceTimes *sTimes = new sequenceTimes[2];
  set_sequenceTimes(sim, sTimes);
  
  /* mcell setup */
  char *tempFile;
  if (d->mcell_file == "") {
    tempFile = mcellTemplate(d->structureFile, d->sheathFile, sim, nProcessors);
  }
  else {
    tempFile = strdup(d->mcell_file.c_str());
  }

  char *logFile;
  if (d->log_file == "") {
    logFile = NULL;
  } else {
    logFile = strdup(d->log_file.c_str());
  }

  mcell_wrap_init(tempFile, nDirs, d->randomSeed + rank, logFile);
  world->msim = (void *)&sim;
  world->is_in_region = is_in_region;

  set_nc_current(&sim, d);

  // Rewrite previous signal files since we later will be appending
  ofstream complexSignalFile(d->complex_signal_filename.c_str(), ios::out);
  complexSignalFile.close();
  //cout << endl << "Writing signal to " << d->signal_filename.c_str() << endl;
  ofstream vertex_signal_file(d->signal_filename.c_str(), ios::out);
  if (vertex_signal_file.fail()) {
    // This really should be checked much earlier! (DR)
    cerr << "Error opening signal file" << d->signal_filename.c_str() << endl;
    exit(-3);
  }
  vertex_signal_file.close();

  //Repeat sequence repeat_number-times
  //TODO Should be able to do repeat time/repeat number with standard acquisition as well
  // probably not currently implemented correctly 
  if (sim.pulse_sequence == PULSE_DWSSFP) {
    sim.endTime = sim.repeat_number*sim.repeat_time;
  }
  if (sim.pulse_sequence == PULSE_MULTI_SPIN_ECHO) {
    sim.endTime += sim.echo_number * sim.echo_spacing;
  }

  int nextUpdate = sim.endTime/100;

  /* start with zero signal */
  int n_species = world->n_species + 1; //0 is for total signal
  complexSignal = new complex<double>[nDirs*n_species];
  vertex_signal = new double[nDirs*n_species];
  gvec_xaxis = new double[3];
  for (int j = 0; j < nDirs*n_species; j++) {
    complexSignal[j] = 0.0;
  }

  /* I don't know what the NSTEPS loop is for, but maybe we'll use it later */
  cout << "\n" << "nDirs: " << nDirs << "\n";
  if (rank == 0) {
    // DR: It seems that some mcell call in the loop below messes up this nice counter.
    cout << "Simulation step completed: ";
    cout.precision(2);
    cout.setf(ios::fixed, ios::floatfield);
    cout.fill(' ');
    cout.width(6);
    //cout << 0.0 << "%";
    cout << endl;
    os.precision(2);
    os.setf(ios::fixed, ios::floatfield);
    os.fill(' ');
    os.width(6);
  }
  for (int istep=0; istep<NSTEPS; istep++) {
    while (currentTime <= sim.endTime - sim.timestep) {
      // Without the following run into issues with running the start of another
      // pulse sequence, probably a better way to fix this but works for now
      if (sim.repeat_number > 1) {
        currentSeqTime = currentTime % sim.repeat_time;
      } else {
        currentSeqTime = currentTime;
      }

      /* diffuse the molecules */
      // Note we need to run initial iteration before anything else to make sure molecules are initialized
      // otherwise things like rf_pulse won't do anything
      mcell_iteration(currentTime/sim.timestep);
#ifdef DEBUG
      count_regions(iteration++);
#endif

      // All of this should probably be in a seperate function
      /* set current gradient */
      switch (sim.pulse_sequence) {
      case PULSE_GRAD_ECHO:
        // Set up magnetization
        if (currentSeqTime == sTimes[0].rampUp1) {
          rf_pulse(90., nDirs, 1, 0);
        }
        if (currentSeqTime < sTimes[0].gOff + sim.gradient_off_time/2) {
          sim.current_gradient_strength =
            CurrentGradientStrength(currentSeqTime, sTimes[0], sim);
        }
        else {
          sim.current_gradient_strength =
            -CurrentGradientStrength(currentSeqTime, sTimes[0], sim);
        }
        break;
      case PULSE_SPIN_ECHO:
        sim.current_gradient_strength =
          CurrentGradientStrength(currentSeqTime, sTimes[0], sim);
        if (currentSeqTime == sTimes[0].rampUp1) {
          rf_pulse(90., nDirs, 1, 0);
        } else if (currentSeqTime >= sTimes[0].gOff + sim.gradient_off_time/2
            && currentSeqTime < sTimes[0].gOff + sim.gradient_off_time/2 + sim.timestep) {
          //rf_180(nDirs);
          rf_pulse(180., nDirs, 0, 0);
        }
        break;
      case PULSE_STEAM_DTI:
        sim.current_gradient_strength =
          CurrentGradientStrength(currentSeqTime, sTimes[0], sim);
        if (currentSeqTime == sTimes[0].rampUp1) {
          rf_pulse(90., nDirs, 1, 0);
        } else if (currentSeqTime >= sTimes[0].gOff + sim.gradient_off_time/2 -
            sim.mixing_time/2 && currentSeqTime < sTimes[0].gOff +
            sim.gradient_off_time/2 - sim.mixing_time/2 + sim.timestep) {
          rf_pulse(90., nDirs, 1, 0);
        } else if (currentSeqTime >= sTimes[0].gOff + sim.gradient_off_time/2 +
            sim.mixing_time/2 && currentSeqTime < sTimes[0].gOff +
            sim.gradient_off_time/2 + sim.mixing_time/2 + sim.timestep) {
          rf_pulse(90., nDirs, 1, 0);
        }
        break;
      case PULSE_DPFG:
        if (d->tensorFileDPFG) {
          float g[2];
          double gvec[2][3];
          double gvec_scaled[2][3];

          // The gradient_direction vector will be scaled appropriately 
          sim.current_gradient_strength = 1.0;
          g[0] = CurrentGradientStrength(currentSeqTime, sTimes[0], sim);
          g[1] = CurrentGradientStrength(currentSeqTime, sTimes[1], sim);
          if (g[0] == 0.0 && g[1] == 0.0) {
            sim.current_gradient_strength = 0.0;
          }
          // Compute actual applied gradient for every pair of gradients, where we assume one
          // gradient is along the x axis. Can be changed later, but there isn't as much value
          // in changing the direction of both gradients
          for (int i = 0; i < gradientDirSize; i++) {
            gvec[0][0] = sim.gx1[i]; gvec[0][1] = sim.gy1[i]; gvec[0][2] = sim.gz1[i];
            gvec[1][0] = sim.gx2[i]; gvec[1][1] = sim.gy2[i]; gvec[1][2] = sim.gz2[i];
            SCALE_VEC3(gvec_scaled[0], g[0], gvec[0]);
            SCALE_VEC3(gvec_scaled[1], g[1], gvec[1]);
            SUM_VEC3(gradient_directions[i], gvec_scaled[0], gvec_scaled[1]);
          }
        } else {
          float g[2];
          double gvec[2][3];
          gvec_xaxis[0] = 1.;
          gvec_xaxis[1] = 0.;
          gvec_xaxis[2] = 0.;

          // The gradient_direction vector will be scaled appropriately 
          sim.current_gradient_strength = 1.0;
          g[0] = CurrentGradientStrength(currentSeqTime, sTimes[0], sim);
          g[1] = CurrentGradientStrength(currentSeqTime, sTimes[1], sim);
          if (g[0] == 0.0 && g[1] == 0.0) {
            sim.current_gradient_strength = 0.0;
          }
          // Compute actual applied gradient for every pair of gradients, where we assume one
          // gradient is along the x axis. Can be changed later, but there isn't as much value
          // in changing the direction of both gradients
          for (int i = 0; i < gradientDirSize; i++) {
            SCALE_VEC3(gvec[0], g[0], gvec_xaxis);
            SCALE_VEC3(gvec[1], g[1], dpfg_gradient_directions[i]);
            SUM_VEC3(gradient_directions[i], gvec[0], gvec[1]);
          }
        }

        if (currentSeqTime == sTimes[0].rampUp1) {
          rf_pulse(90., nDirs, 1, 0);
        } else if ((currentSeqTime >= sTimes[0].gOff + sim.gradient_off_time/2
             && currentSeqTime < sTimes[0].gOff + sim.gradient_off_time/2 + sim.timestep)
            || (currentSeqTime >= sTimes[1].gOff + sim.gradient_off_time/2
                && currentSeqTime < sTimes[1].gOff + sim.gradient_off_time/2 + sim.timestep)) {
          //rf_180(nDirs);
          rf_pulse(180., nDirs, 0, 0);
        }
        break;
      case PULSE_DWSSFP:
        sim.current_gradient_strength = CurrentGradientStrength(currentSeqTime, sTimes[0], sim);
        if (currentSeqTime == sTimes[0].rampUp1) {
          // Random precession angles from unbalanced gradients or slice selection
          //rf_pulse(0., nDirs, 2., 1);
          SaveSignalFile(sim, d, nDirs, currentTime);
          rf_pulse(sim.flip_angle, nDirs, 1, 0);
          SaveSignalFile(sim, d, nDirs, currentTime);
        }
        break;
      case PULSE_MULTI_SPIN_ECHO:
        sim.current_gradient_strength = CurrentGradientStrength(currentSeqTime, sTimes[0], sim);
        if (currentSeqTime == 0) {
          rf_pulse(90., nDirs, 1, 0);
        } else if (currentSeqTime >= sTimes[0].gOff + sim.gradient_off_time/2
            && currentSeqTime < sTimes[0].gOff + sim.gradient_off_time/2 + sim.timestep) {
          rf_pulse(180., nDirs, 0, 0);
        }
        //Only check for rf pulse if more than 1 echo, otherwise need to handle math errors when echo_time = 0
        if (sim.echo_number > 1) {
          if ((currentSeqTime > sim.echo_time) && ((currentSeqTime - sim.echo_time)%sim.echo_spacing >= sim.echo_spacing/2) && ((currentSeqTime - sim.echo_time)%sim.echo_spacing < sim.echo_spacing/2 + sim.timestep)) {
            rf_pulse(180, nDirs, 0, 0);
          }
        }
        break;

      }

      if (rank == 0 && currentTime > nextUpdate) {
        os.str("");
        os.width(6);
        os << (int)(100000.0*currentTime/(1.0*sim.endTime))/1000.0 << "%";
        t_now = time(NULL);
        t_end = t_now+(t_now - t_start)/(double)currentTime*
          (double)(sim.endTime-currentTime);
        diffusionStatistics(sim, maxDx, maxDy, maxDz);
        os << " Max Diffusion (x,y,z): ("
           << maxDx << "," << maxDy << "," << maxDz << ")";
        os << " ERT: " << (long)difftime(t_end, t_start)
             << " ETC: " << (long)difftime(t_end, t_now) << " ";
        for (auto && p : os.str()) cout << "\b\b";
        cout << os.str();
        cout.flush();
        nextUpdate += sim.endTime/100;
      }

      update_nc_current(&sim, currentSeqTime);
      mcell_update_spins(&sim, nDirs, gradient_directions, rank);
        
      if (sim.pulse_sequence == PULSE_DWSSFP) {
        if (currentSeqTime == sTimes[0].gOff) {
          //Update magnetization to reflect precession accrued during diffusion gradient
          update_magnetization(sim, nDirs);
        }
        //If last time step of a DWSSFP sequence, apply all relaxation
        //TODO Check if this is actually the right time step, maybe dt not subtracted?
        if (currentSeqTime == sim.repeat_time-sim.timestep) {
          update_relaxation(sim, currentSeqTime, nDirs);
          mcell_reset_origpos(&sim, nDirs);
        }
      } else {
        // Should figure out how to update relaxation less frequently (not every time step)
        // similar to how we do DWSSFP. Can this be tied directly to rf_pulses and be ok?
        if (sim.pulse_sequence == PULSE_DPFG) {
          if (currentSeqTime == sTimes[0].gOff || currentSeqTime == sTimes[0].end || currentSeqTime == sTimes[1].gOff || currentSeqTime == sTimes[1].end) {
            update_magnetization(sim, nDirs);
          }
        }
        else if (currentSeqTime == sTimes[0].gOff || currentSeqTime == sTimes[0].end) {
          //Update magnetization to reflect precession accrued during diffusion gradient
          update_magnetization(sim, nDirs);
        }
#if 0
        if (sim.use_relaxation) {
          update_relaxation(sim, sim.timestep, nDirs);
        }
#endif
      }

      if (sim.pulse_sequence == PULSE_MULTI_SPIN_ECHO) {
        if ((currentSeqTime == sim.echo_time) || ((currentSeqTime > sim.echo_time) && (currentSeqTime-sim.echo_time)%sim.echo_spacing == 0))  {
          if(SaveSignalFile(sim, d, nDirs, currentTime)) {
            cerr << " Error saving signal file " << endl;
            exit(-1);
          }
        }
      }
      
      currentTime += sim.timestep;
      B0.t = currentTime;

    } /* end while (currentTime < sim.endTime) loop */

    if (sim.pulse_sequence != PULSE_DWSSFP) {
      if (sim.use_relaxation) {
        update_relaxation(sim, B0.t-B0.tp, nDirs);
      }
    }

    if (rank == 0) {
      for (auto && p : os.str()) cout << "\b";
      cout << "\b\b\b\b\b\b\b";
      cout.width(6);
      cout << 100.0 << "%" << endl;
      cout.flush();
    }

    diffusionStatistics(sim, currentTime);

    /* compute signal now */
    // TODO Need to make sure magnetization is in sync with phase_shift at this point. Is it ok
    // or maybe necessary to update_magnetization here? Seems ok for now, only case where its a problem
    // is if something weird happens to phase_shift I think.
    update_magnetization(sim, nDirs);
    switch (sim.pulse_sequence) {
    case PULSE_DWSSFP:
    case PULSE_MULTI_SPIN_ECHO:
      //DWSSFP and MULTI_SPIN_ECHO have different signal output
      break;
    case PULSE_GRAD_ECHO:
    case PULSE_SPIN_ECHO:
    case PULSE_DPFG:
    case PULSE_STEAM_DTI:
      if(SaveSignalFile(sim, d, nDirs, currentTime)) {
        cerr << " Error saving signal file " << endl;
        exit(-1);
      }
    }

    mcell_wrap_close();
    if (tempFile) {
      if (rank != 0) {
        remove(tempFile);
      }
    }
  } // steps loop
}

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
    case 'h': printUsage(cerr); exit(0);
      break;
    default:
      cerr << "Unknown option: " << (char)next_option << endl;
      printUsage(cerr);
      exit(1);
    }
  }
}

int main(int argc, char *argv[]) {
  progname = argv[0];
  parseOptions(argc, argv);
  if (optind==argc) { printUsage(cerr); exit(0); }
  input_data d;
  parse_input(&d, argv[1]);
  simulation(&d);
  return 0;
}
