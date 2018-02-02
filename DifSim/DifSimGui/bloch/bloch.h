/*<BLOCH.h>***************************************************************/

#ifndef _BLOCH_HEADER_
#define _BLOCH_HEADER_

/* header file for bloch */

/* Generic includes */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/fcntl.h> 

#ifndef MATCH_STR
#  define MATCH_STR(s1,s2) ((strcmp(s1,s2))?(0):(1))
#endif

#define MAXINT 32752.0
#define GAM    26747.5    /* gamma in 1/(s-G) */
#define BIGRF  100000
#define INT_TO_RAD (M_PI/(MAXINT)) /* Integer to radian conversion */
#define BIGVAL 1.e+10

#define RF_INVERT  0
#define RF_EXCITE  1
#define RF_REFOCUS 2

/* my includes */

#include "alloc.h"
#include "coords.h"
#include "vector.h"
#include "matrix.h"
#include "ioutils.h"
#include "excite.h"

/* rf includes */
#include "vsechgp.h"
#include "genmb.h"

#endif /* _BLOCH_HEADER_ */
