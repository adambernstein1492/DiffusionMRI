/*<ioutils.h>***************************************************************/

#ifndef _LRF_IO_HEADER_
#define _LRF_IO_HEADER_

#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <sys/fcntl.h> 
#include <fcntl.h> 

#include <sys/time.h>
#include <sys/times.h>

#ifndef YES
#define YES 1 
#endif
#ifndef NO
#define NO 0 
#endif
#ifndef MAYBE
#define MAYBE -1
#endif

#define TMP_EXIT \
{fprintf(stderr,"Temporary exit (from routine %s) ...\n",program_name); exit(0);}

/********** Warnings *********************/

#ifndef FATAL_ERROR
#define FATAL_ERROR(msg) \
        {fprintf(stderr,"%s\a\n in file: %s at line %d\n ... exiting\n", \
                 msg,__FILE__,__LINE__);exit(1);}
#endif

#ifndef NEAR_FATAL_ERROR
#define NEAR_FATAL_ERROR(msg) \
        {fprintf(stderr,"%s\a\n in file: %s at line %d\n ... returning\n", \
                 msg,__FILE__,__LINE__);return(NULLVAL);}
#endif

#ifndef NEAR_FATAL_ZERROR
#define NEAR_FATAL_ZERROR(msg) \
        {fprintf(stderr,"%s\a\n in file: %s at line %d\n ... returning\n", \
                 msg,__FILE__,__LINE__);return(ZNULLVAL);}
#endif

#ifndef WARNING
#define WARNING(msg) \
        {fprintf(stderr,"%s\a\n in file: %s at line %d\n", \
                 msg,__FILE__,__LINE__);}
#endif

#ifndef MATCH_STR
#  define MATCH_STR(s1,s2) ((strcmp((s1),(s2)))?(0):(1))
#endif
#ifndef MATCH_STR_N
#  define MATCH_STR_N(s1,s2,n) ((strncmp((s1),(s2),(n)))?(0):(1))
#endif
#ifndef MATCH_STR_STR
#  define MATCH_STR_STR(s1,s2) ((strstr((s1),(s2))==NULL)?(0):(1))
#endif

/*== File macros ========================================================*/

#define RESET_FILE(fp) fseek((fp),0,SEEK_SET)
#define SET_FILE_LOC(fp,loc) fseek((fp),(loc),SEEK_SET)
#define PRINT_FILE_LOC(fp) printf("Current file location = %i\n",ftell((fp)))

#define FILE_LINES(file,nlines)                          \
  do{ int len=512;                                       \
      char *cbuf,nc;                                     \
      cbuf = (char *)malloc((unsigned) len * sizeof(char)); \
      RESET_FILE((file));\
      nlines = 0; \
      while ( fgets(cbuf,len,(file)) != NULL ) { \
	(nlines) ++ ; \
      } ; \
  } while(0)

/* Set file "fp" location to "line" */

#define SET_FILE(fp,line) \
 do{ char str[MAXPROC] ; \
     int  nlines = 0; \
     if((fp)==NULL) { \
       printf("SET_FILE: File not open!\n"); \
     } else {\
       RESET_FILE((fp)); \
       while ((fgets(str,MAXPROC,(fp)) != NULL) && ((nlines++)<((line)-2)) ); \
       printf("SET_FILE: line = %i, nlines = %i\n",(line),nlines); \
       fseek((fp),ftell((fp)),SEEK_SET) ; \
     }\
 } while(0)

/*******************************/

#ifndef BETWEEN
#define BETWEEN(x,a,b)  (((x) >= (a)) && ((x)<=(b)))
#endif

#define SET_LIMS(lims,xmin,xmax)  ((lims)[0]=(xmin),(lims)[1]=(xmax))

static char QandA[256] ;	/* question and answer string, 
				   for use with sinput */

void cox_sleep( int msec );
int file_exists( char *fname );
int file_lines( char *fname );

double pinput(char *pname, double pdef, double *plims) ;
char *sinput(char *sname, char *sdef, char *slims[], int ns);
bool logqyn(char *question, char *qdef) ;

void  read_iqm(long int *ier,char fname[],long int *size,char arr[]) ;
void write_iqm(long int *ier,char fname[],long int *size,char arr[]) ;
void printRange(double *x,int nx, char *xname) ;

#endif
