#ifndef MCELL_UTIL
#define MCELL_UTIL

#include <stdio.h>

#define COUNT_OF(arr) (sizeof((arr))/sizeof((arr[0])))

/**********************************************************************
* Definitions for the infinite array whose size may grow. *
*
***********************************************************************/

#define BLOCK_SIZE 10000

/** struct infinite_double_array
    Used to hold information for an infinite array of doubles.
*/
struct infinite_double_array{
	/* the data  for this block */
	double data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_double_array *next;
};

/** struct infinite_int_array
    Used to hold information for an infinite array of integers.
*/
struct infinite_int_array{
	/* the data  for this block */
	int data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_int_array *next;
};

/** struct infinite_uint_array
    Used to hold information for an infinite array of unsigned integers.
*/
struct infinite_uint_array{
	/* the data  for this block */
	unsigned int data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_uint_array *next;
};

/** struct infinite_longlong_array
    Used to hold information for an infinite array of long long integers.
*/
struct infinite_longlong_array{
	/* the data  for this block */
	long long data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_longlong_array *next;
};

/** struct infinite_string_array
    Used to hold information for an infinite array of strings.
*/
struct infinite_string_array{
	/* the data  for this block */
	char *data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_string_array *next;
};

/** struct infinite_pointer_array
    Used to hold information for an infinite array of pointers.
*/
struct infinite_pointer_array{
	/* the data  for this block */
        /* array of pointers */
	void *data[BLOCK_SIZE];

	/* pointer to the next array. */
	struct infinite_pointer_array *next;
};

/* Initializes the infinite array. 
   Initialization of the field "data" should be performed
   separately to the value corresponding to the data type of the field array*/
#define ia_init(array_ptr)	{(array_ptr)->next = NULL;}

struct iteration_counter
{
    long long *iterations; /* array of iteration numbers, should be 
                              memory allocated */
    int max_iterations; /* size of the above array */
    int n_iterations;   /* number of the filled positions in an above array */
};

struct string_buffer
{
    char **strings;   /* array of strings, should be memory allocated */
    int max_strings;  /* size of the above array */
    int n_strings;    /* number of the filled positions in an above array */
};

double ia_double_get(struct infinite_double_array *array_ptr, int idx);
void ia_double_store(struct infinite_double_array *array_ptr, int idx, double data_to_store);
int ia_int_get(struct infinite_int_array *array_ptr, int idx);
void ia_int_store(struct infinite_int_array *array_ptr, int idx, int data_to_store);
unsigned int ia_uint_get(struct infinite_uint_array *array_ptr, int idx);
void ia_uint_store(struct infinite_uint_array *array_ptr, int idx, unsigned int data_to_store);
long long ia_longlong_get(struct infinite_longlong_array *array_ptr, long long idx);
void ia_longlong_store(struct infinite_longlong_array *array_ptr, long long idx, long long data_to_store);
char* ia_string_get(struct infinite_string_array *array_ptr, int idx);
void ia_string_store(struct infinite_string_array *array_ptr, int idx, char *data_to_store);
void *ia_pointer_get(struct infinite_pointer_array *array_ptr, int idx);
void ia_pointer_store(struct infinite_pointer_array *array_ptr, int idx, void *data_to_store);


struct bit_array
{
  int nbits;
  int nints;
  /* Bit array data runs off the end of this struct */
};

struct bit_array* new_bit_array(int bits);
struct bit_array* duplicate_bit_array(struct bit_array *old);
int get_bit(struct bit_array* ba, int idx);
void set_bit(struct bit_array *ba, int idx, int value);
void set_bit_range(struct bit_array *ba,int idx1,int idx2,int value);
void set_all_bits(struct bit_array *ba,int value);
void bit_operation(struct bit_array *ba,struct bit_array *bb,char op);
int count_bits(struct bit_array *ba);
void print_bit_array(FILE *F, struct bit_array *ba);
void free_bit_array(struct bit_array *ba);


int bisect(double *list,int n,double val);
int bisect_near(double *list,int n,double val);
int bisect_high(double *list,int n,double val);
int bin(double *list,int n,double val);


int distinguishable(double a,double b,double eps);

int is_abbrev(char *abbrev,char *full);
int is_reverse_abbrev(char *abbrev,char *full);

struct void_list
{
  struct void_list *next;
  void *data;
};

struct void_list* void_list_sort(struct void_list *vl);
struct void_list* void_list_sort_by(struct void_list *vl,int (*leq)(void*,void*));
void remove_one_duplicate(struct void_list *sorted);
int remove_both_duplicates(struct void_list **head);
void delete_void_list(struct void_list *head);


int void_array_search(void **array,int n,void *to_find);
int void_ptr_compare(void const *v1, void const *v2);

unsigned int *allocate_uint_array(int size, unsigned int value);
void **allocate_ptr_array(int size);
void free_ptr_array(void **pa, int count);

struct num_expr_list;
void free_num_expr_list(struct num_expr_list *nel);
void uniq_num_expr_list(struct num_expr_list *nel);

int is_dir(char const *path);
int is_writable_dir(char const *path);
int mkdirs(char const *path);
int make_parent_dir(char const *path);

FILE *open_file(char const *fname, char const *mode);
int get_basename(char const *filepath, char **basename);
int get_dirname(char const *filepath, char **dirname);

#if ! (defined(__ICC) || defined(__INTEL_COMPILER))
double erfcinv(double v);
#define erfinv(x) erfcinv(1-(x))
#endif

int poisson_dist(double lambda,double p);

void byte_swap(void *data, int size);

/* This function analyzes the string and checks
   whether the string contains wildcards (*, ?,[,]).
   Returns 1 if wildcard is found, and 0 - otherwise). */
int contain_wildcard(char *teststring);

int feral_strlenn(char *feral,int n);
int is_feral_nabbrev(char *feral,int n,char *tame);
char* feral_strstrn(char *tame_haystack,char *feral_needle,int n);
int is_wildcard_match(char *wild,char *tame);

int dir_exists(char const *filename);
int initialize_iteration_counter(struct iteration_counter *cntr, int max_iters);
int destroy_iteration_counter(struct iteration_counter *cntr);
int add_to_iteration_counter_monotonic(struct iteration_counter *cntr, long long iter);
int add_to_iteration_counter(struct iteration_counter *cntr, long long iter);
int initialize_string_buffer(struct string_buffer *sb, int maxstr);
int destroy_string_buffer(struct string_buffer *sb);
int add_string_to_buffer(struct string_buffer *sb, char *str);

/*******************************************************************
 Pointer hashes

   Pointer hashes were originally written for the table lookups in
   macromolecules -- primarily for the cooperative rates.  Essentially,
   they are pointer -> pointer hash tables.  There is no restriction on
   the type of pointer used for the key, but a hash value must be
   supplied along with the pointer whenever performing any operation
   that requires a key.  Neither the key, nor the value pointer are
   ever dereferenced or freed by the pointer hash code.

      Usage:

        struct pointer_hash hash;
        if (pointer_hash_init(&hash)) { fail }
        pointer_hash_add(&hash, key1, hash1, value1);
        pointer_hash_add(&hash, key2, hash2, value2);
        ..
        pointer_hash_add(&hash, keyN, hashN, valueN);

        ..

        void *value1 = pointer_hash_lookup(&hash, key1, hash1);

        ..

        pointer_hash_destroy(&hash);

  Note: The pointer hash allocates some memory, which will be orphaned
  if pointer_hash_destroy is not called on the hash before it goes out
  of scope.

  The hash table collision strategy implemented by this structure is a
  streaming strategy, rather than a chaining strategy.  This means that
  when we try to insert a value into the n-th bin, if the n-th bin is
  already occupied, we will move to the n+1-th bin (modulo table size)
  until we find an insertion point.  When we do a lookup, therefore, we
  need to scan forward until we either find our key, or until we find
  an empty bin.  We can ensure that at least one of these conditions is
  met by always keeping the table size greater than the number of
  entries.

  Currently, the pointer hash is tuned to keep the table size between 2
  and 4 times the number of items it contains.
*******************************************************************/
struct pointer_hash
{
  int           num_items;    /* num items in table */
  int           table_size;   /* size of table */
  unsigned int *hashes;       /* hash values for each entry */
  void const  **keys;         /* keys for each entry */
  void        **values;       /* values for each entry */
};

/* Initialize a pointer hash to a given initial size.  Returns 0 on
 * success. */
int pointer_hash_init(struct pointer_hash *ht,
                      int size);

/* Quickly clear all values from a pointer hash.  Does not free any
 * memory. */
void pointer_hash_clear(struct pointer_hash *ht);

/* Destroy a pointer hash, freeing all memory associated with it. */
void pointer_hash_destroy(struct pointer_hash *ht);

/* Manually resize a pointer hash to have at least 'new_size' bins.
 * New size may exceed requested size. */
int pointer_hash_resize(struct pointer_hash *ht,
                        int new_size);

/* Add a value to a pointer hash.  If a previous item was added for
 * that key, the new value will replace the old value. */
int pointer_hash_add(struct pointer_hash *ht,
                     void const *key,
                     unsigned int keyhash,
                     void *value);

/* Look up a value in a pointer hash.  Returns NULL if no item was
 * found, or if the value associated with the key was NULL. */
#define pointer_hash_lookup(ht, key, keyhash) \
            pointer_hash_lookup_ext(ht, key, keyhash, NULL)

/* Look up a value in a pointer hash.  Returns NULL if no item was
 * found, or if the value associated with the key was NULL. */
void *pointer_hash_lookup_ext(struct pointer_hash const *ht,
                              void const *key,
                              unsigned int keyhash,
                              void *default_value);

/* Remove a value from a pointer hash.  Returns 0 if the item was
 * successfully removed, or 0 if the item was not found.
 */
int pointer_hash_remove(struct pointer_hash *ht,
                        void const *key,
                        unsigned int keyhash);

/*******************************************************************
 Inline min/max functions
*******************************************************************/

static inline double min2d(double x, double y)
{
  return (x < y) ? x : y;
}

static inline int min2i(int x, int y)
{
  return (x < y) ? x : y;
}

static inline double max2d(double x, double y)
{
  return (x > y) ? x : y;
}

static inline int max2i(int x, int y)
{
  return (x > y) ? x : y;
}

static inline double min3d(double x, double y, double z)
{
    return (z<y) ? ((z < x)?z:x) : ((y<x)?y:x);
}

static inline int min3i(int x, int y, int z)
{
    return (z<y) ? ((z < x)?z:x) : ((y<x)?y:x);
}

static inline double max3d(double x, double y, double z)
{
    return (z>y) ? ((z > x)?z:x) : ((y>x)?y:x);
}

static inline int max3i(int x, int y, int z)
{
    return (z>y) ? ((z > x)?z:x) : ((y>x)?y:x);
}

static inline long long min3ll(long long x,long long y,long long z)
{
    return (z<y) ? ((z < x)?z:x) : ((y<x)?y:x);
    
}

static inline long long max3ll(long long x,long long y,long long z)
{
    return (z>y) ? ((z > x)?z:x) : ((y>x)?y:x);
}    


/* Return minimum value from the array of N doubles */ 
static inline double minNd(double *array, int N)
{
    double smallest;
    N-=2;
    for (smallest = array[N+1]; N >= 0; N--)
    {
      if (array[N] < smallest) smallest=array[N];
    }
    return smallest;
}

/* Return minimum value from the array of N integers */ 
static inline int minNi(int *array, int N)
{
    int smallest;
    N-=2;
    for (smallest = array[N+1]; N >= 0; N--)
    {
      if (array[N] < smallest) smallest=array[N];
    }
    return smallest;
}

#endif
