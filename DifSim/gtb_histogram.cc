#include "gtb_histogram.h"

#include <stdio.h>

gtb_histogram *gtb_histogram_alloc(size_t n) {
  gtb_histogram *h = (gtb_histogram *) malloc(sizeof(gtb_histogram));
  h->size = n;
  h->bin = (int *) calloc(n, sizeof(int));
  h->range = (double *) calloc((n+1), sizeof(double));
  return h;
}

void gtb_histogram_free(gtb_histogram *h) {
  free(h->bin);
  free(h->range);
  free(h);
}

int gtb_histogram_set_uniform_ranges(gtb_histogram *h, double min, double max) {
  if (!(h->size > 0) || !h->bin || !h->range) {
    /* something isn't initialized properly */
    return 1;
  }
  for (int i = 0; i < h->size; i++) {
    h->range[i] = i*((max - min)/h->size) + min;
  }
  h->range[h->size] = max;
  return 0;
}

int gtb_histogram_increment(gtb_histogram *h, double x) {
  if (!(h->size > 0) || !h->bin || !h->range) {
    /* something isn't initialized properly */
    return 1;
  }
  for (int i = 0; i < h->size; i++) {
    if (x >= h->range[i] && x < h->range[i+1]) {
      h->bin[i]++;
      return 0;
    }
  }
  return 1;
}

double gtb_histogram_get_bin_mean(gtb_histogram *h, size_t i) {
  if (i < 0 || i >= h->size) {
    return 0.0;
  }
  return (h->range[i] + h->range[i+1])/2.0;
}
  
int gtb_histogram_get_value(gtb_histogram *h, size_t i) {
  if (i < 0 || i >= h->size) {
    return -1;
  }
  return h->bin[i];
}

const char *gtb_histogram_print(gtb_histogram *h) {
  size_t s = 24*h->size*sizeof(char);
  char *r = (char *)malloc(s);
  char *r_tmp = r;
  int w = 0;
  int w_tmp;
  for (int i = 0; i < h->size; i++) {
    w_tmp = snprintf(r_tmp, s-w, "%3.7f %10d\n", (h->range[i] + h->range[i+1])/2.0, h->bin[i]);
    r_tmp += w_tmp;
    w += w_tmp;
  }
  return (const char *)r;
}
