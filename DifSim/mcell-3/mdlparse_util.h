#ifndef MDLPARSE_UTIL_H
#define MDLPARSE_UTIL_H

#include "vector.h"
#include "mcell_structs.h"
#include "mdlparse_aux.h"

#ifndef __GNUC__
#ifndef __attribute__
#define __attribute__(x) /* empty */
#endif
#endif

#if __GNUC__ < 3
#ifndef __attribute__
#define __attribute__(x) /* empty */
#endif
#endif

/* ====================================
 * Caveat lector: Most of the functions in this file (and the corresponding .c
 * file) are designed to be called from the grammar, and take care of cleaning
 * up memory allocation.  If you are in doubt about whether you can still use
 * any pointers you've passed in check the function source code to be certain.
 * In most cases, the answer is "no".
 * ====================================
 */

/* Strips the quotes from a quoted string, freeing the original. */
char *mdl_strip_quotes(struct mdlparse_vars *mpvp, char *in);

/* Concatenates two strings, freeing the original strings. */
char *mdl_strcat(struct mdlparse_vars *mpvp, char *s1, char *s2);

/* Duplicates a string. */
char *mdl_strdup(struct mdlparse_vars *mpvp, char const *s1);

/* Display a warning message about something encountered during the
 * parse process. */
void mdl_warning(struct mdlparse_vars *mpvp, char const *fmt, ...)
    __attribute__((format (printf, 2, 3)));

/* Find an include file based on the path of the currently parsed file
 */
char *mdl_find_include_file(struct mdlparse_vars  *mpvp,
                            char const *path,
                            char const *cur_path);

/* Check that the speficied file mode string is valid for an fopen statement.
 */
int mdl_valid_file_mode(struct mdlparse_vars *mpvp,
                        char *mode);

/* Error-reporting version of log */
int mdl_expr_log(struct mdlparse_vars *mpvp, double in, double *out);

/* Error-reporting version of log10 */
int mdl_expr_log10(struct mdlparse_vars *mpvp, double in, double *out);

/* Error-reporting version of MOD */
int mdl_expr_mod(struct mdlparse_vars *mpvp,
                 double in,
                 double divisor,
                 double *out);

/* Error-reporting version of division operator */
int mdl_expr_div(struct mdlparse_vars *mpvp,
                 double in,
                 double divisor,
                 double *out);

/* Error-reporting version of exponentiation operator */
int mdl_expr_pow(struct mdlparse_vars *mpvp,
                 double in,
                 double exponent,
                 double *out);

/* Get a uniform random number */
double mdl_expr_rng_uniform(struct mdlparse_vars *mpvp);

/* Round a value off to n digits */
double mdl_expr_roundoff(struct mdlparse_vars *mpvp,
                         double in,
                         int ndigits);

/* Turn a string into a double */
int mdl_expr_string_to_double(struct mdlparse_vars *mpvp,
                              char *str,
                              double *out);

/* Add a new file handle to the symbol table. */
struct sym_table *mdl_new_filehandle(struct mdlparse_vars *mpvp,
                                     char *name);

/* Process an fopen statement, opening a new file handle. */
int mdl_fopen(struct mdlparse_vars *mpvp,
              struct sym_table *filesym,
              char *name,
              char *mode);

/* Process an fclose statement, closing an existing file handle. */
int mdl_fclose(struct mdlparse_vars *mpvp,
               struct sym_table *filesym);

/* Create a new double argument for a printf argument list. */
struct arg *mdl_new_printf_arg_double(struct mdlparse_vars *mpvp, double d);

/* Create a new string argument for a printf argument list. */
struct arg *mdl_new_printf_arg_string(struct mdlparse_vars *mpvp,
                                      char const *str);

/* Expands C-style escape sequences in the string. */
char *mdl_expand_string_escapes(struct mdlparse_vars *mpvp, char const *in);

/* printf-like formatting of MDL arguments.  Prints to the defined err_file. */
int mdl_printf(struct mdlparse_vars *mpvp, char *fmt, struct arg *arg_head);

/* fprintf-like formatting of MDL arguments. */
int mdl_fprintf(struct mdlparse_vars *mpvp,
                struct file_stream *filep,
                char *fmt,
                struct arg *arg_head);

/* (expression-form) sprintf-like formatting of MDL arguments. */
char *mdl_string_format(struct mdlparse_vars *mpvp,
                        char *fmt,
                        struct arg *arg_head);

/* sprintf-like formatting of MDL arguments. */
int mdl_sprintf(struct mdlparse_vars *mpvp,
                struct sym_table *assign_var,
                char *fmt,
                struct arg *arg_head);

/* strtime-like formatting of current time. */
int mdl_fprint_time(struct mdlparse_vars *mpvp,
                    struct sym_table *filep,
                    char *fmt);

/* strtime-like formatting of current time.  Prints to the defined err_file. */
void mdl_print_time(struct mdlparse_vars *mpvp, char *fmt);

/* Generate a num_expr_list containing the numeric values from start to end,
 * incrementing by step. */
int mdl_generate_range(struct mdlparse_vars *mpvp,
                       struct num_expr_list_head *list,
                       double start,
                       double end,
                       double step);

/* Generate a numeric list containing a single value. */
int mdl_generate_range_singleton(struct mdlparse_vars *mpvp,
                                 struct num_expr_list_head *lh,
                                 double value);

/* Add a value to a numeric list. */
int mdl_add_range_value(struct mdlparse_vars *mpvp,
                        struct num_expr_list_head *lh,
                        double value);

#ifdef DEBUG
/* Display a human-readable representation of the specified array. */
void mdl_debug_dump_array(struct num_expr_list *nel);
#else
#define mdl_debug_dump_array(x) /* do nothing */
#endif

/* Create a 3-D vector from a numeric array. */
struct vector3 *mdl_point(struct mdlparse_vars *mpvp,
                          struct num_expr_list_head *vals);

/* Create a 3-D vector equal to s*[1, 1, 1] for some scalar s. */
struct vector3 *mdl_point_scalar(struct mdlparse_vars *mpvp,
                                 double val);

/* Get a named variable if it exists, or create it if it doesn't. */
struct sym_table *mdl_get_or_create_variable(struct mdlparse_vars *mpvp,
                                             char *name);

/* Assign a "double" value to a variable, freeing any previous value. */
int mdl_assign_variable_double(struct mdlparse_vars *mpvp,
                               struct sym_table *sym,
                               double value);

/* Assign a string value to a variable, freeing any previous value. */
int mdl_assign_variable_string(struct mdlparse_vars *mpvp,
                               struct sym_table *sym,
                               char const *value);

/* Assign an array value to a variable, freeing any previous value. */
int mdl_assign_variable_array(struct mdlparse_vars *mpvp,
                              struct sym_table *sym,
                              struct num_expr_list *value);

/* Assign one variable value to another variable, freeing any previous value.
 */
int mdl_assign_variable(struct mdlparse_vars *mpvp,
                        struct sym_table *sym,
                        struct sym_table *value);

/* Set all notification levels to a particular value. */
void mdl_set_all_notifications(struct volume *vol, byte notify_value);

/* Set the iteration report frequency. */
int mdl_set_iteration_report_freq(struct mdlparse_vars *mpvp, long long interval);

/* Set all warning levels to a particular value. */
void mdl_set_all_warnings(struct volume *vol, byte warning_level);

/* Set the lifetime warning threshold. */
int mdl_set_lifetime_warning_threshold(struct mdlparse_vars *mpvp, long long lifetime);

/* Set the missed reaction warning threshold. */
int mdl_set_missed_reaction_warning_threshold(struct mdlparse_vars *mpvp, double rxfrac);

/* Set the global timestep for the simulation. */
int mdl_set_time_step(struct mdlparse_vars *mpvp, double step);

/* Set the maximum timestep for the simulation. */
int mdl_set_max_time_step(struct mdlparse_vars *mpvp, double step);

/* Set the global space step for the simulation. */
int mdl_set_space_step(struct mdlparse_vars *mpvp, double step);

/* Set the number of iterations for the simulation. */
int mdl_set_num_iterations(struct mdlparse_vars *mpvp, long long numiters);

/* Set the number of radial directions. */
int mdl_set_num_radial_directions(struct mdlparse_vars *mpvp, int numdirs);

/* Set the number of radial subdivisions. */
int mdl_set_num_radial_subdivisions(struct mdlparse_vars *mpvp, int numdivs);

/* Set the effector grid density. */
int mdl_set_grid_density(struct mdlparse_vars *mpvp, double density);

/* Set the number of times to place any particular macromolecule. */
int mdl_set_complex_placement_attempts(struct mdlparse_vars *mpvp,
                                       double attempts);

/* Schedule an asynchronous checkpoint. */
int mdl_set_realtime_checkpoint(struct mdlparse_vars *mpvp,
                                long duration,
                                int cont_after_cp);

/* Set the input checkpoint file to use. */
int mdl_set_checkpoint_infile(struct mdlparse_vars *mpvp,
                              char *name);

/* Set the output checkpoint file to use. */
int mdl_set_checkpoint_outfile(struct mdlparse_vars *mpvp,
                               char *name);

/* Set the number of iterations between checkpoints. */
int mdl_set_checkpoint_interval(struct mdlparse_vars *mpvp,
                                long long iters);

/* Set the partitioning in a particular dimension. */
int mdl_set_partition(struct mdlparse_vars *mpvp,
                      int dim,
                      struct num_expr_list_head *head);

/* Starts creating a new object.  This function has side effects, and must be
 * paired with a call to mdl_finish_object(mdlpvp) */
struct sym_table *mdl_start_object(struct mdlparse_vars *mpvp,
                                   char *name);

/* "Finishes" a new object, undoing the state changes that occurred when the
 * object was "started".
 */
void mdl_finish_object(struct mdlparse_vars *mpvp);

/* Adds the first element to an empty object list. */
void mdl_object_list_singleton(struct mdlparse_vars *mpvp,
                               struct object_list *head,
                               struct object *objp);

/* Adds an element to an object list. */
void mdl_add_object_to_list(struct mdlparse_vars *mpvp,
                            struct object_list *head,
                            struct object *objp);

/* Find an existing object or print an error message if the object isn't found.
 */
struct sym_table *mdl_existing_object(struct mdlparse_vars *mpvp,
                                      char *name);

/* Find a list of existing objects matching a particular wildcard. */
struct sym_table_list *mdl_existing_objects_wildcard(struct mdlparse_vars *mpvp,
                                                     char *wildcard);

/* Find an existing region or print an error message if it isn't found. */
struct sym_table *mdl_existing_region(struct mdlparse_vars *mpvp,
                                      struct sym_table *obj_symp,
                                      char *name);

/* Find an existing molecule species, or print an error message if it isn't
 * found. */
struct sym_table *mdl_existing_molecule(struct mdlparse_vars *mpvp,
                                        char *name);

/* Turn a single symbol into a singleton symbol list. */
struct sym_table_list *mdl_singleton_symbol_list(struct mdlparse_vars *mpvp,
                                                 struct sym_table *sym);

/* Find an existing molecule species, and return it in a singleton list, or
 * print an error message if it isn't found. */
struct sym_table_list *mdl_existing_molecule_list(struct mdlparse_vars *mpvp,
                                                  char *name);

/* Find a list of all molecule species matching the specified wildcard.  Print
 * an error message if it doesn't match any. */
struct sym_table_list *mdl_existing_molecules_wildcard(struct mdlparse_vars *mpvp,
                                                       char *wildcard);

/* Find an existing macromolecule species, or print an error message if it
 * isn't found, or isn't a macromolecule. */
struct sym_table *mdl_existing_macromolecule(struct mdlparse_vars *mpvp,
                                             char *name);

/* Find an existing surface molecule species, or print an error message if it
 * isn't found, or isn't a surface molecule. */
struct sym_table *mdl_existing_surface_molecule(struct mdlparse_vars *mpvp,
                                                char *name);

/* Find an existing surface class species, or print an error message if it
 * isn't found, or isn't a surface class. */
struct sym_table *mdl_existing_surface_class(struct mdlparse_vars *mpvp,
                                             char *name);

/* Find a named variable if it exists, or print an error if it does not. */
struct sym_table *mdl_existing_variable(struct mdlparse_vars *mpvp,
                                        char *name);

/* Find an existing array symbol, or print an error message if it isn't found.
 */
struct sym_table *mdl_existing_array(struct mdlparse_vars *mpvp, char *name);

/* Find a named numeric variable if it exists, or print an error message if it
 * isn't found. */
struct sym_table *mdl_existing_double(struct mdlparse_vars *mpvp, char *name);

/* Find a named string variable if it exists, or print an error message if it
 * isn't found. */
struct sym_table *mdl_existing_string(struct mdlparse_vars *mpvp, char *name);

/* Get a named numeric or array variable if it exists.  Print an error message
 * if it isn't found. */
struct sym_table *mdl_existing_num_or_array(struct mdlparse_vars *mpvp,
                                            char *name);

/* Find an existing named reaction pathway or molecule, or print an error
 * message if it isn't found. */
struct sym_table *mdl_existing_rxn_pathname_or_molecule(struct mdlparse_vars *mpvp,
                                                        char *name);

/* Find an existing reaction pathway or release pattern, or print an error
 * message if it isn't found, or if the name could refer to either a release
 * pattern or a reaction pathway.
 */
struct sym_table *mdl_existing_release_pattern_or_rxn_pathname(struct mdlparse_vars *mpvp,
                                                              char *name);

/* Find an existing molecule or object, or print an error message if it isn't
 * found, or if the name could refer to either type of object. */
struct sym_table *mdl_existing_molecule_or_object(struct mdlparse_vars *mpvp, char *name);

/* Find an existing file stream or print an error message if it isn't found. */
struct sym_table *mdl_existing_file_stream(struct mdlparse_vars *mpvp,
                                           char *name);

/* Find all mesh objects (polygons and boxes) that match a given wildcard. */
struct sym_table_list *mdl_meshes_by_wildcard(struct mdlparse_vars *mpvp,
                                              char *wildcard);

/* Apply a translation to the given transformation matrix. */
void mdl_transform_translate(struct mdlparse_vars *mpvp,
                             double (*mat)[4],
                             struct vector3 *xlat);

/* Apply a scale to the given transformation matrix. */
void mdl_transform_scale(struct mdlparse_vars *mpvp,
                         double (*mat)[4],
                         struct vector3 *scale);

/* Apply a rotation to the given transformation matrix. */
int mdl_transform_rotate(struct mdlparse_vars *mpvp,
                         double (*mat)[4],
                         struct vector3 *axis,
                         double angle);

/* Deep copy an object. */
int mdl_deep_copy_object(struct mdlparse_vars *mpvp,
                         struct object *dst_obj,
                         struct object *src_obj);

/* Prepare a region description for use in the simulation by creating a
 * membership bitmask on the region object. */
int mdl_normalize_elements(struct mdlparse_vars *mpvp, struct region *reg,int existing);

/* Finalizes the polygonal structure of the box, normalizing all regions. */
int mdl_triangulate_box_object(struct mdlparse_vars *mpvp,
                               struct sym_table *box_sym,
                               struct polygon_object *pop,
                               double box_aspect_ratio);

/* Clean up the regions on an object, eliminating any removed walls. */
void mdl_remove_gaps_from_regions(struct object *ob);

/* Check that the specified diffusion constant is valid, correcting it if
 * appropriate. */
int mdl_check_diffusion_constant(struct mdlparse_vars *mpvp, double *d);

/* Finish the creation of a molecule, undoing any state changes we made during
 * the creation of the molecule. */
void mdl_finish_molecule(struct mdlparse_vars *mpvp,
                         struct species *mol);

/* Finish the creation of a series of molecules, undoing any state changes we
 * made during the creation of the molecules. */
void mdl_finish_molecules(struct mdlparse_vars *mpvp,
                          struct species_list_item *mols);

/* Populate a species list with a single species. */
int mdl_species_list_singleton(struct mdlparse_vars *mpvp,
                               struct species_list *list,
                               struct species *spec);

/* Add a single species to a species list. */
int mdl_add_to_species_list(struct mdlparse_vars *mpvp,
                            struct species_list *list,
                            struct species *spec);

/* Start parsing the innards of a release site. */
int mdl_start_release_site(struct mdlparse_vars *mpvp,
                           struct sym_table *symp,
                           int shape);

/* Finish parsing the innards of a release site. */
struct object *mdl_finish_release_site(struct mdlparse_vars *mpvp,
                                       struct sym_table *symp);

/* Set the location of a release site. */
void mdl_set_release_site_location(struct mdlparse_vars *mpvp,
                                   struct release_site_obj *rsop,
                                   struct vector3 *location);

/* Validate a release site. */
int mdl_is_release_site_valid(struct mdlparse_vars *mpvp, struct release_site_obj *rsop);

/* Set the geometry for a particular release site to be a region expression. */
int mdl_set_release_site_geometry_region(struct mdlparse_vars *mpvp,
                                         struct release_site_obj *rsop,
                                         struct object *objp,
                                         struct release_evaluator *re);

/* Set the geometry for a particular release site to be an entire object. */
int mdl_set_release_site_geometry_object(struct mdlparse_vars *mpvp,
                                         struct release_site_obj *rsop,
                                         struct object *objp);

/* Create a new "release on region" expression term. */
struct release_evaluator *mdl_new_release_region_expr_term(struct mdlparse_vars *mpvp,
                                                           struct sym_table *my_sym);

/* Set the geometry for a particular release site to be a region expression. */
struct release_evaluator *mdl_new_release_region_expr_binary(struct mdlparse_vars *mpvp,
                                                             struct release_evaluator *reL,
                                                             struct release_evaluator *reR,
                                                             int op);

/* Set the molecule to be released from this release site. */
int mdl_set_release_site_molecule(struct mdlparse_vars *mpvp,
                                  struct release_site_obj *rsop,
                                  struct species_opt_orient *mol_type);

/* Set the diameter of a release site. */
int mdl_set_release_site_diameter(struct mdlparse_vars *mpvp,
                                  struct release_site_obj *rsop,
                                  double diam);

/* Set the diameter of the release site along the X, Y, and Z axes. */
int mdl_set_release_site_diameter_array(struct mdlparse_vars *mpvp,
                                        struct release_site_obj *rsop,
                                        int n_diams,
                                        struct num_expr_list *diams,
                                        double factor);

/* Set the diameters of the release site along the X, Y, and Z axes from a
 * variable, either scalar or vector. */
int mdl_set_release_site_diameter_var(struct mdlparse_vars *mpvp,
                                      struct release_site_obj *rsop,
                                      double factor,
                                      struct sym_table *symp);

/* Set the release probability for a release site. */
int mdl_set_release_site_probability(struct mdlparse_vars *mpvp,
                                     struct release_site_obj *rsop,
                                     double prob);

/* Set the release pattern to be used by a particular release site. */
int mdl_set_release_site_pattern(struct mdlparse_vars *mpvp,
                                 struct release_site_obj *rsop,
                                 struct sym_table *pattern);

/* Set the molecule positions for a LIST release. */
int mdl_set_release_site_molecule_positions(struct mdlparse_vars *mpvp,
                                            struct release_site_obj *rsop,
                                            struct release_single_molecule_list *list);

/* Create a mew single molecule release position for a LIST release site. */
struct release_single_molecule *mdl_new_release_single_molecule(struct mdlparse_vars *mpvp,
                                                                struct species_opt_orient *mol_type,
                                                                struct vector3 *pos);

/* Populates a list with a single LIST release molecule descriptor. */
void mdl_release_single_molecule_singleton(struct mdlparse_vars *mpvp,
                                           struct release_single_molecule_list *list,
                                           struct release_single_molecule *mol);

/* Adds a release molecule descriptor to a list. */
void mdl_add_release_single_molecule_to_list(struct mdlparse_vars *mpvp,
                                             struct release_single_molecule_list *list,
                                             struct release_single_molecule *mol);

/* Set a constant release quantity from this release site, in units of
 * molecules. */
void mdl_set_release_site_constant_number(struct mdlparse_vars *mpvp,
                                          struct release_site_obj *rsop,
                                          double num);

/* Set a gaussian-distributed release quantity from this release site, in units
 * of molecules. */
void mdl_set_release_site_gaussian_number(struct mdlparse_vars *mpvp,
                                          struct release_site_obj *rsop,
                                          double mean,
                                          double stdev);

/* Set a release quantity from this release site based on a fixed concentration
 * in a sphere of a gaussian-distributed diameter with a particular mean and
 * std. deviation. */
void mdl_set_release_site_volume_dependent_number(struct mdlparse_vars *mpvp,
                                                  struct release_site_obj *rsop,
                                                  double mean,
                                                  double stdev,
                                                  double conc);

/* Set a release quantity from this release site based on a fixed concentration
 * within the release-site's area. */
int mdl_set_release_site_concentration(struct mdlparse_vars *mpvp,
                                       struct release_site_obj *rsop,
                                       double conc);

/* Set a release quantity from this release site based on a fixed density
 * within the release-site's area. */
int mdl_set_release_site_density(struct mdlparse_vars *mpvp,
                                 struct release_site_obj *rsop,
                                 double dens);

/* Set an item to be the sole element of a vertex list. */
void mdl_vertex_list_singleton(struct mdlparse_vars *mpvp,
                               struct vertex_list_head *head,
                               struct vertex_list *item);

/* Append a vertex to a list. */
void mdl_add_vertex_to_list(struct mdlparse_vars *mpvp,
                            struct vertex_list_head *head,
                            struct vertex_list *item);

/* Allocate an item for a vertex list. */
struct vertex_list *mdl_new_vertex_list_item(struct mdlparse_vars *mpvp,
                                             struct vector3 *vertex);

/* Set an item to be the sole element of an element connection list. */
void mdl_element_connection_list_singleton(struct mdlparse_vars *mpvp,
                                           struct element_connection_list_head *head,
                                           struct element_connection_list *item);

/* Append an element connection to a list. */
void mdl_add_element_connection_to_list(struct mdlparse_vars *mpvp,
                                        struct element_connection_list_head *head,
                                        struct element_connection_list *item);

/* Create an element connection (essentially a triplet of vertex indices). */
struct element_connection_list *mdl_new_element_connection(struct mdlparse_vars *mpvp,
                                                           struct num_expr_list_head *indices);

/* Create a tetrahedral element connection (essentially a quadruplet of vertex
 * indices). */
struct element_connection_list *mdl_new_tet_element_connection(struct mdlparse_vars *mpvp,
                                                               struct num_expr_list_head *indices);

/* Create a new polygon list object. */
struct polygon_object *mdl_new_polygon_list(struct mdlparse_vars *mpvp,
                                            struct sym_table *sym,
                                            int n_vertices,
                                            struct vertex_list *vertices,
                                            int n_connections,
                                            struct element_connection_list *connections);

/* Finalize the polygon list, cleaning up any state updates that were made when
 * we started creating the polygon. */
int mdl_finish_polygon_list(struct mdlparse_vars *mpvp,
                            struct sym_table *symp);

/* Check a box or polygon list object for degeneracy. */
int mdl_check_degenerate_polygon_list(struct mdlparse_vars *mpvp,
                                      struct object *objp);

/* Create a new voxel list object. */
struct voxel_object *mdl_new_voxel_list(struct mdlparse_vars *mpvp,
                                        struct sym_table *sym,
                                        int n_vertices,
                                        struct vertex_list *vertices,
                                        int n_connections,
                                        struct element_connection_list *connections);

/* Create a new box object, with particular corners. */
struct polygon_object *mdl_new_box_object(struct mdlparse_vars *mpvp,
                                          struct sym_table *sym,
                                          struct vector3 *llf,
                                          struct vector3 *urb);

/* Finalize the box object, cleaning up any state updates that were made when
 * we started creating the box. */
int mdl_finish_box_object(struct mdlparse_vars *mpvp,
                          struct sym_table *symp);

/* Create a named region on an object. */
struct region *mdl_create_region(struct mdlparse_vars *mpvp,
                                 struct object *objp,
                                 char *name);

/* Get a region on an object, creating it if it does not exist yet. */
struct region *mdl_get_region(struct mdlparse_vars *mpvp,
                              struct object *objp,
                              char *name);

/* Begin construction of a region on an existing object. */
int mdl_start_existing_obj_region_def(struct mdlparse_vars *mpvp,
                                      struct sym_table *obj_symp);

/* Append elements to an element list. */
void mdl_add_elements_to_list(struct mdlparse_vars *mpvp,
                              struct element_list_head *list,
                              struct element_list *head,
                              struct element_list *tail);

/* Marks elements as being excluded, rather than included (the default). */
void mdl_set_elements_to_exclude(struct mdlparse_vars *mpvp,
                                 struct element_list *els);

/* Create a new element list for a region description. */
struct element_list *mdl_new_element_list(struct mdlparse_vars *mpvp,
                                          unsigned int begin,
                                          unsigned int end);

/* Create a new element list for a region description based on a side name. */
struct element_list *mdl_new_element_side(struct mdlparse_vars *mpvp,
                                          unsigned int side);

/* Create a new element list for a "previous region" include/exclude statement.
 */
struct element_list *mdl_new_element_previous_region(struct mdlparse_vars *mpvp,
                                                     struct object *objp,
                                                     struct region *rp_container,
                                                     char *name_region_referent,
                                                     int exclude);

/* Allocate a new region element list item for an include/exclude PATCH
 * statement. */
struct element_list *mdl_new_element_patch(struct mdlparse_vars *mpvp,
                                           struct polygon_object *poly,
                                           struct vector3 *llf,
                                           struct vector3 *urb,
                                           int exclude);

/* Set the elements for a region, normalizing the region if it's on a polygon
 * list object. */
int mdl_set_region_elements(struct mdlparse_vars *mpvp,
                            struct region *rgn,
                            struct element_list *elements,
                            int normalize_now);

/* Create a new named reaction pathway name structure. */
struct sym_table *mdl_new_rxn_pathname(struct mdlparse_vars *mpvp,
                                       char *name);

/* Adds children to a meta-object, aggregating counts of walls and vertices
 * from the children into the specified parent.  The children should already
 * have their parent pointers set. */
void mdl_add_child_objects(struct mdlparse_vars *mpvp,
                           struct object *parent,
                           struct object *child_head,
                           struct object *child_tail);

/* Adds an effector (or list of effectors) to a region.  These effectors will
 * be placed on the surface at initialization time. */
void mdl_add_effector_to_region(struct mdlparse_vars *mpvp,
                                struct region *rgn,
                                struct eff_dat_list *lst);

/* Set the surface class of this region, possibly inheriting the viz_value.  */
void mdl_set_region_surface_class(struct mdlparse_vars *mpvp,
                                  struct region *rgn,
                                  struct sym_table *scsymp);

/* Set the VIZ_VALUE for this region. */
void mdl_set_region_region_viz_value(struct mdlparse_vars *mpvp,
                                     struct region *rgn,
                                     int viz_value);

/****************************************************************
 * Reaction output
 ***************************************************************/

/* Finalizes a reaction data output block, checking for errors, and allocating
 * the output buffer. */
int mdl_output_block_finalize(struct mdlparse_vars *mpvp,
                              struct output_block *obp);

/* Create a new output set for reaction output. */
struct output_set* mdl_new_output_set(struct mdlparse_vars *mpvp,
                                      char *comment,
                                      int exact_time);

/* Populate an output set. */
struct output_set *mdl_populate_output_set(struct mdlparse_vars *mpvp,
                                           struct output_set *os,
                                           struct output_column *col_head,
                                           int file_flags,
                                           char *outfile_name);

/* Construct and add an output block to the world. */
int mdl_add_reaction_output_block_to_world(struct mdlparse_vars *mpvp,
                                           int buffer_size,
                                           struct output_times_inlist *otimes,
                                           struct output_set_list *osets);

/* Joins two subtrees into a reaction data output expression tree, with a
 * specified operation. */
struct output_expression* mdl_join_oexpr_tree(struct mdlparse_vars *mpvp,
                                              struct output_expression *left,
                                              struct output_expression *right,
                                              char oper);

/* Converts an output expression tree generated from a wildcard into a
 * summation expression tree. */
struct output_expression *mdl_sum_oexpr(struct mdlparse_vars *mpvp,
                                        struct output_expression *expr);

/* Creates a constant output expression for reaction data output. */
struct output_expression *mdl_new_oexpr_constant(struct mdlparse_vars *mpvp,
                                                 double value);

/* Generates a reaction data output expression from the first count syntax form
 * (simple molecule, unquoted, no orientation). */
struct output_expression *mdl_count_syntax_1(struct mdlparse_vars *mpvp,
                                             struct sym_table *what,
                                             struct sym_table *where,
                                             int hit_spec,
                                             int count_flags);

/* Generates a reaction data output expression from the second count syntax
 * form (simple molecule, unquoted, orientation in braces). */
struct output_expression *mdl_count_syntax_2(struct mdlparse_vars *mpvp,
                                             struct sym_table *mol_type,
                                             short orient,
                                             struct sym_table *where,
                                             int hit_spec,
                                             int count_flags);

/* Generates a reaction data output expression from the third count syntax form
 * (quoted string, possibly a wildcard, possibly an oriented molecule). */
struct output_expression *mdl_count_syntax_3(struct mdlparse_vars *mpvp,
                                             char *what,
                                             struct sym_table *where,
                                             int hit_spec,
                                             int count_flags);

/* Generate a reaction data output expression from the macromolecule "subunit"
 * syntax variant. */
struct output_expression *mdl_count_syntax_macromol_subunit(struct mdlparse_vars *mpvp,
                                                            struct complex_species *macromol,
                                                            struct species_opt_orient *master_orientation,
                                                            struct species_opt_orient *subunit,
                                                            struct macro_relation_state *relation_states,
                                                            struct sym_table *location);

/* Prepare a single count expression for inclusion in an output set. */
int mdl_single_count_expr(struct mdlparse_vars *mpvp,
                          struct output_column_list *list,
                          struct output_expression *expr,
                          char *custom_header);

/****************************************************************
 * Viz output
 ***************************************************************/

/* Build a new VIZ output block, containing parameters for an output set for
 * visualization. */
int mdl_new_viz_output_block(struct mdlparse_vars *mpvp);

/* Finalize a new VIZ output block, ensuring all required parameters were set,
 * and doing any postprocessing necessary for runtime. */
int mdl_finish_viz_output_block(struct mdlparse_vars *mpvp,
                                struct viz_output_block *vizblk);

/* Require the mode of the specified VIZ output block to be pre-MCell 3 (i.e.
 * neither of the DREAMM_V3 modes.)  It's impossible to create a valid DREAMM
 * output using the old-style notation, so this makes the impossibility more
 * explicit for user-friendliness. */
int mdl_require_old_style_viz(struct mdlparse_vars *mpvp, int mode);

/* Set the mode for a new VIZ output block. */
int mdl_set_viz_mode(struct mdlparse_vars *mpvp,
                     struct viz_output_block *vizblk,
                     int mode);

/* Set the mesh format for a new VIZ output block. */
int mdl_set_viz_mesh_format(struct mdlparse_vars *mpvp,
                            struct viz_output_block *vizblk,
                            int format);

/* Set the molecule format for a new VIZ output block. */
int mdl_set_viz_molecule_format(struct mdlparse_vars *mpvp,
                                struct viz_output_block *vizblk,
                                int format);

/* Set the filename prefix for a new VIZ output block. */
int mdl_set_viz_filename_prefix(struct mdlparse_vars *mpvp,
                                struct viz_output_block *vizblk,
                                char *filename);

/* Set the molecule filename prefix for an old-style VIZ output block. */
int mdl_set_viz_molecule_filename_prefix(struct mdlparse_vars *mpvp,
                                         struct viz_output_block *vizblk,
                                         char *filename);

/* Set the object filename prefix for an old-style VIZ output block. */
int mdl_set_viz_object_filename_prefix(struct mdlparse_vars *mpvp,
                                       struct viz_output_block *vizblk,
                                       struct sym_table *sym,
                                       char *filename);

/* Error-checking wrapper for a specified visualization state. */
int mdl_viz_state(struct mdlparse_vars *mpvp,
                  int *target,
                  double value);

/* Sets a flag on all of the listed objects, requesting that they be
 * visualized. */
int mdl_set_viz_include_meshes(struct mdlparse_vars *mpvp,
                               struct viz_output_block *vizblk,
                               struct sym_table_list *list,
                               int viz_state);

/* Sets the viz state of a particular mesh object, indicating whether it should
 * be visualized. */
int mdl_set_viz_include_mesh_state(struct mdlparse_vars *mpvp,
                                   struct viz_output_block *vizblk,
                                   struct sym_table *obj,
                                   int viz_state);

/* Sets a flag on a viz block, requesting that all meshes be visualized. */
int mdl_set_viz_include_all_meshes(struct mdlparse_vars *mpvp,
                                   struct viz_output_block *vizblk,
                                   int viz_state);

/* Sets a flag on all of the listed species, requesting that they be visualized. */
int mdl_set_viz_include_molecules(struct mdlparse_vars *mpvp,
                                  struct viz_output_block *vizblk,
                                  struct sym_table_list *list,
                                  int viz_state);

/* Sets a flag on a viz block, requesting that all species be visualized. */
int mdl_set_viz_include_all_molecules(struct mdlparse_vars *mpvp,
                                      struct viz_output_block *vizblk,
                                      int viz_state);

/* Adds some new mesh output frames to a list. */
int mdl_new_viz_mesh_frames(struct mdlparse_vars *mpvp,
                            struct viz_output_block *vizblk,
                            struct frame_data_list_head *frames,
                            int time_type,
                            int mesh_item_type,
                            struct num_expr_list_head *timelist);

/* Adds some new molecule output frames to a list. */
int mdl_new_viz_mol_frames(struct mdlparse_vars *mpvp,
                           struct viz_output_block *vizblk,
                           struct frame_data_list_head *frames,
                           int time_type,
                           int mol_item_type,
                           struct num_expr_list_head *timelist);

/* Adds some new output frames to a list. */
int mdl_new_viz_frames(struct mdlparse_vars *mpvp,
                       struct viz_output_block *vizblk,
                       struct frame_data_list_head *frames,
                       int time_type,
                       int type,
                       struct num_expr_list_head *times);

/* Build a list of times for VIZ output, one timepoint per iteration in the
 * simulation. */
int mdl_new_viz_all_times(struct mdlparse_vars *mpvp,
                          struct num_expr_list_head *list);

/* Build a list of iterations for VIZ output, one for each iteration in the
 * simulation. */
int mdl_new_viz_all_iterations(struct mdlparse_vars *mpvp,
                               struct num_expr_list_head *list);

/* Set the viz_state value for an object and all of its children. */
int mdl_set_object_viz_state_by_name(struct mdlparse_vars *mpvp,
                                     struct viz_output_block *vizblk,
                                     struct sym_table *obj_symp,
                                     int viz_state);

/* Set the viz_state value for a molecular species. */
int mdl_set_molecule_viz_state(struct mdlparse_vars *mpvp,
                               struct viz_output_block *vizblk,
                               struct species *specp,
                               int viz_state);

/* Set the viz_state for a particular region. */
int mdl_set_region_viz_state(struct mdlparse_vars *mpvp,
                             struct viz_output_block *vizblk,
                             struct region *rp,
                             int viz_state);

/* Allocate a block of data for Rex's custom visualization mode. */
struct rk_mode_data *mdl_new_rk_mode_var(struct mdlparse_vars *mpvp,
                                         struct num_expr_list_head *values,
                                         struct vector3 *direction);

/****************************************************************
 * Volume output
 ***************************************************************/

/* Create a new volume output request. */
struct volume_output_item *mdl_new_volume_output_item(struct mdlparse_vars *mpvp,
                                                      char *filename_prefix,
                                                      struct species_list *molecules,
                                                      struct vector3 *location,
                                                      struct vector3 *voxel_size,
                                                      struct vector3 *voxel_count,
                                                      struct output_times *ot);

/* Create new default output timing for volume output. */
struct output_times *mdl_new_output_times_default(struct mdlparse_vars *mpvp);

/* Create new "step" output timing for volume output. */
struct output_times *mdl_new_output_times_step(struct mdlparse_vars *mpvp,
                                               double step);

/* Create new "iteration list" output timing for volume output. */
struct output_times *mdl_new_output_times_iterations(struct mdlparse_vars *mpvp,
                                                     struct num_expr_list_head *iters);

/* Create new "time list" output timing for volume output. */
struct output_times *mdl_new_output_times_time(struct mdlparse_vars *mpvp,
                                               struct num_expr_list_head *times);

/****************************************************************
 * Release patterns
 ***************************************************************/

/* Create a new release pattern.  There must not yet be a release pattern with
 * the given name. */
struct sym_table *mdl_new_release_pattern(struct mdlparse_vars *mpvp,
                                          char *name);

/* Fill in the details of a release pattern. */
int mdl_set_release_pattern(struct mdlparse_vars *mpvp,
                            struct sym_table *rpat_sym,
                            struct release_pattern *rpat_data);

/****************************************************************
 * Molecules
 ***************************************************************/

/* Is the given name a valid macromolecule name? */
int mdl_valid_complex_name(struct mdlparse_vars *mpvp, char *name);

/* Create a new species.  There must not yet be a molecule or named reaction
 * pathway with the supplied name. */
struct sym_table *mdl_new_molecule(struct mdlparse_vars *mpvp, char *name);

/* Assemble a molecule species from its component pieces. */
struct species *mdl_assemble_mol_species(struct mdlparse_vars *mpvp,
                                         struct sym_table *sym,
                                         double D_ref,
                                         double D,
                                         int is_2d,
                                         double time_step,
                                         int target_only,
                                         double max_step_length,
                                         double relax_t1,
                                         double relax_t2,
                                         double lambda_1,
                                         double lambda_2,
                                         double lambda_3);

/****************************************************************
 * Reactions, surface classes
 ***************************************************************/

/* Check whether the reaction rate is valid. */
int mdl_valid_rate(struct mdlparse_vars *mpvp,
                   struct reaction_rate *rate);

/* Set the reactant/product list to contain a single item. */
int mdl_reaction_player_singleton(struct mdlparse_vars *mpvp,
                                  struct species_opt_orient_list *list,
                                  struct species_opt_orient *spec);

/* Add a single item to a reactant/product player list. */
int mdl_add_reaction_player(struct mdlparse_vars *mpvp,
                            struct species_opt_orient_list *list,
                            struct species_opt_orient *spec);

/* Set a reaction rate from a variable. */
int mdl_reaction_rate_from_var(struct mdlparse_vars *mpvp,
                               struct reaction_rate *rate,
                               struct sym_table *symp);

/* Set a complex reaction rate. */
int mdl_reaction_rate_complex(struct mdlparse_vars *mpvp,
                              struct reaction_rate *rate,
                              struct sym_table *symp,
                              char *tbl);

/* Assemble a standard reaction from its component parts. */
struct rxn *mdl_assemble_reaction(struct mdlparse_vars *mpvp,
                                  struct species_opt_orient *reactants,
                                  struct species_opt_orient *surface_class,
                                  struct reaction_arrow *react_arrow,
                                  struct species_opt_orient *products,
                                  struct reaction_rates *rate,
                                  struct sym_table *pathname);

/* Assemble a surface reaction from its component parts. */
struct rxn *mdl_assemble_surface_reaction(struct mdlparse_vars *mpvp,
                                          int reaction_type,
                                          struct species *surface_class,
                                          struct sym_table *reactant_sym,
                                          short orient);

/* Assemble a concentration clamp reaction from its component parts. */
struct rxn *mdl_assemble_concentration_clamp_reaction(struct mdlparse_vars *mpvp,
                                                      struct species *surface_class,
                                                      struct sym_table *mol_sym,
                                                      short orient,
                                                      double conc);

/* Start a surface class declaration. */
void mdl_start_surface_class(struct mdlparse_vars *mpvp,
                             struct sym_table *symp);

/* Finish a surface class declaration.  Undoes side effects from
 * mdl_start_surface_class. */
void mdl_finish_surface_class(struct mdlparse_vars *mpvp,
                              struct sym_table *symp);

/* Create a new effector data for surface molecule initialization. */
struct eff_dat *mdl_new_effector_data(struct mdlparse_vars *mpvp,
                                      struct species_opt_orient *eff_info,
                                      double quant);


/****************************************************************
 * Macromolecules
 ***************************************************************/


/* Allocate a new relation state structure for a rule table (presently, either
 * a rate table, or a subunit counting table). */
struct macro_relation_state *mdl_assemble_complex_relation_state(struct mdlparse_vars *mpvp,
                                                                 int rel_idx,
                                                                 int invert,
                                                                 struct species_opt_orient *mol);

/* Assemble a subunit assignment for one or more subunits within a complex. */
struct macro_subunit_assignment *mdl_assemble_complex_subunit_assignment(struct mdlparse_vars *mpvp,
                                                                         struct macro_subunit_spec *su,
                                                                         struct species_opt_orient *spec);

/* Assemble a complex geometry structure. */
struct macro_geometry *mdl_assemble_complex_geometry(struct mdlparse_vars *mpvp,
                                                     struct macro_topology *topo,
                                                     struct num_expr_list_head *coords,
                                                     struct vector3 *pos);

/* Checks the parsed relationships for errors, freeing them if errors are
 * found, and passing them on, otherwise.  rels are freed if they are invalid.
 */
int mdl_validate_complex_relationships(struct mdlparse_vars *mpvp,
                                       struct macro_topology *topo,
                                       struct macro_relationship *rels);

/* Assemble a complex relationship. */
struct macro_relationship *mdl_assemble_complex_relationship(struct mdlparse_vars *mpvp,
                                                             struct macro_topology *topo,
                                                             char *name,
                                                             struct num_expr_list_head *rel);

/* Validate the rate rules, freeing them if invalid, and passing them on if
 * valid.
 */
int mdl_validate_complex_rates(struct mdlparse_vars *mpvp,
                               struct macro_rate_ruleset *rates);

/* Assemble a macromolecular complex rate rule set. */
struct macro_rate_ruleset *mdl_assemble_complex_ruleset(struct mdlparse_vars *mpvp,
                                                        char *name,
                                                        struct macro_rate_rule *rules);

/* Assemble a macromolecular complex rate rule. */
struct macro_rate_rule *mdl_assemble_complex_rate_rule(struct mdlparse_vars *mpvp,
                                                       struct macro_rate_clause *clauses,
                                                       double rate);

/* Assemble a macromolecular rate rule clause. */
struct macro_rate_clause *mdl_assemble_complex_rate_rule_clause(struct mdlparse_vars *mpvp,
                                                                struct macro_relationship *rels,
                                                                char *relation_name,
                                                                int invert,
                                                                struct species_opt_orient *target);

/* Validate the geometry of a macromolecular complex, freeing it if it is
 * invalid. */
int mdl_validate_complex_geometry(struct mdlparse_vars *mpvp,
                                  struct macro_topology *topo,
                                  struct macro_geometry *geom);

/* Allocate a new topology data structure for a macromolecule. */
struct macro_topology *mdl_assemble_topology(struct mdlparse_vars *mpvp,
                                             struct num_expr_list_head *dims);

/* Allocate and populate a new component in a subunit coordinate specification.
 */
struct macro_subunit_spec *mdl_assemble_subunit_spec_component(struct mdlparse_vars *mpvp,
                                                               int from,
                                                               int to);

/* Assemble a complex species, adding it to the symbol table.  */
int mdl_assemble_complex_species(struct mdlparse_vars *mpvp,
                                 char *name,
                                 struct macro_topology *topo,
                                 struct macro_subunit_assignment *assignments,
                                 struct macro_geometry *geom,
                                 struct macro_relationship *rels,
                                 struct macro_rate_ruleset *rates);

/* Postprocess the parsed reactions, moving them to the reaction hash table,
 * and transferring information from the pathway structures to a more compact,
 * runtime-optimized form.
 */
int prepare_reactions(struct mdlparse_vars *mpvp);

void free_vertex_list(struct vertex_list *vlp);
#endif
