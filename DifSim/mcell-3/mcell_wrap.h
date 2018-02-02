#ifndef MCELL_WRAP_H
#define MCELL_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif

  void mcell_iteration(long long it_time);
  int mcell_wrap_close();
  void mcell_wrap_init(char *file_name, int nDirs, int seed, char *log_file);
  struct subvolume* find_subvolume_wrap(struct vector3 *loc);
  struct subvolume* next_subvol_wrap(struct vector3 *start,struct vector3 *end,struct subvolume *sv);
  int collide_wall_wrap(struct vector3 *start,struct vector3 *move,struct wall *face, double *t,struct vector3 *hitpt,int update_move);

  double wrap_rng();
  int bisect_wrap(double *list,int n,double val);
  double collide_sv_time_wrap(struct vector3 *here,struct vector3 *move,struct subvolume *sv);

#ifdef __cplusplus
}
#endif

#endif /* MCELL_WRAP_H */
