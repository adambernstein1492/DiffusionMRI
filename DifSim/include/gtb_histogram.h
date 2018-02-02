#ifndef GTB_HISTOGRAM_H
#define GTB_HISTOGRAM_H

#include <stdlib.h>

typedef struct gtb_histogram_struct {
  size_t size;       /* for holding the number of bins */
  int *bin;          /* for holding counts of each range */
  double *range;     /* for holding the range boundaries */
} gtb_histogram;


gtb_histogram *gtb_histogram_alloc(size_t n);
void gtb_histogram_free(gtb_histogram *h);
int gtb_histogram_set_uniform_ranges(gtb_histogram *h, double min, double max);
int gtb_histogram_increment(gtb_histogram *h, double x);
double gtb_histogram_get_bin_mean(gtb_histogram *h, size_t i);
int gtb_histogram_get_value(gtb_histogram *h, size_t i);
const char *gtb_histogram_print(gtb_histogram *h);
 
#endif /* GTB_HISTOGRAM_H */
