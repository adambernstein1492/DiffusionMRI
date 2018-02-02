#include "ioutils.h"

void printRange(double *x,int nx, char *xname)
{
  int i;

  double xmax = -1.e10;
  double xmin =  1.e10;
  for (i=0; i<nx; i++) {
    if (x[i]<xmin) xmin = x[i];
    if (x[i]>xmax) xmax = x[i];
  }

  fprintf(stderr,"min(%s) = %lf, max(%s) = %lf\n",xname,xmin,xname,xmax);
  return ;
}

/*-----------------------------------------------------------------
   Sleep a given # of milliseconds (uses the Unix select routine)
------------------------------------------------------------------*/

void cox_sleep( int msec )
{
   struct timeval tv ;
   if( msec <= 0 ) return ;
   tv.tv_sec  = msec/1000 ;
   tv.tv_usec = (msec%1000)*1000 ;
   select( 1 , NULL,NULL,NULL , &tv ) ;
   return ;
}

/*===<file_exists.c>==========================================================

     Title: FILE EXISTS

   Purpose: Test if a file exists 

      Call: int file_exists(char *fname)
*/
int file_exists( char *fname )
{
  FILE *outfile;
        
  outfile = fopen (fname,"r");
  if (outfile == NULL)
    return (0);
  else 
    fclose (outfile);
  return (1);
}

/*===<file_lines.c>============================================================

     Title: FILE LINES determination

   Purpose: Determine the number of lines in a file

      Call: int file_lines ( char *fname )
*/
int file_lines( char *fname )
{
  FILE *fin ;
  int  nlines=0 ;
  int  len = 512 ;
  char *cbuf,nc ;

  if( (fin = fopen( fname , "r" )) == NULL ){
    printf("Couldn't open header file %s\n" , fname ) ;
    return (-1) ;
  }

  cbuf = (char *)malloc((unsigned) len * sizeof(char));
  while ( fgets(cbuf,len,fin) != NULL ) {
    nlines ++ ;
  }

  fclose( fin ) ;

  return ( nlines ) ;
}

/*===<pinput.c>========================================================

     Title: Parameter INPUT

   Purpose: Prompt user for numerical input

      Call: double pinput(char *pname, double pdef)
      
     Input: char *pname = parameter name for prompt
            double pdef  = parameter default value

    Return: double value of input

   History: 02_11_15 - lrf

*/
double pinput(char *pname, double pdef, double *plims)
{
  char  line[100];
  int   nn;
  double dval;

  int mistakes = 0;
  int too_many_mistakes = 10;

  do {

    printf(" Enter %s [%g] {%g,%g}: ",pname,pdef,plims[0],plims[1]);
    fgets(line,sizeof(line),stdin);
    nn = sscanf(line,"%lf",&dval);

    if(nn<0) dval = pdef;

    if(!BETWEEN(dval,plims[0],plims[1])) {
      printf("\n ***** Parameter '%s' must be between %lf and %lf\a\n"
	     "(current value = %lf\n)",
	     pname,plims[0],plims[1],dval) ;
      mistakes++ ;
      if(mistakes>too_many_mistakes) {
	printf("You've had your chance!  Later ...\a\n");
	exit(0);
      }
    }

  }  while(!BETWEEN(dval,plims[0],plims[1])) ;

  return(dval) ;
}

/*===<logqyn.c>========================================================

     Title: LOGical Question: Yes or No

   Purpose: Prompt user for yes/no answer

      Call: bool logqyn(char *question, char *qdef)
      
     Input: char *question = query to respond to

    Return: YES or NO

   History: 02_12_30 - lrf

*/
bool logqyn(char *question, char *qdef)
{
  char line[100], sval[100] ;
  int  nn, answer ;

  do {

    printf(" %s? [%s]: ",question,qdef);
    fgets(line,sizeof(line),stdin);
    nn = sscanf(line,"%s",sval);

    if(nn<0) strcpy(sval,qdef);

         if( !strncmp(sval,"y",1) || !strncmp(sval,"Y",1) ) answer = YES   ;
    else if( !strncmp(sval,"n",1) || !strncmp(sval,"N",1) ) answer = NO    ;
    else                                                    answer = MAYBE ;

  } while (answer!=YES && answer!=NO) ;

  return(answer) ;
}

/*===<sinput.c>========================================================

     Title: String INPUT

   Purpose: Prompt user for string input

      Call: char sinput(char *sval, char *sname, char *sdef)
      
     Input: char *sname = string name for prompt
            char sdef   = string default value
            char slims  = allowable values to enter

    Return: char value of input

   History: 02_11_15 - lrf

     Notes: Need to put in allowable values
*/
char *sinput(char *sname, char *sdef, char *slims[], int ns)
{
  char i,line[100];
  char *sval;
  bool match = 0;

  bool debug_local = 0;

  if(debug_local){
    int i;
    fprintf(stderr,"sinput: \n");
    for (i=0; i<ns; i++) 
      fprintf(stderr,"\tslims[%i] = %s\n",i,slims[i]);
  }

  do {
    fprintf(stderr," Enter %s {%s}: ",sname,sdef);
    fgets(line,sizeof(line),stdin);
    if(strlen(line)==1) {
      sval = sdef;
      match = 1;
    } else {
      for (i=0; i<ns; i++) {
	if (!strncmp(line,slims[i],1)) {
	  sval = line;
	  match = 1;
	}
      }
    }
  } while (match==0) ;

  return(sval) ;
}

/**<read_iqm.c>***********************************************
 
    Title: READ I/Q channels
 
  Purpose: Read complex data

     Call: void read_iqm(ier,fname,size,arr)
             long int *ier;
	     long int *size;
	     char fname[],arr[];

   ier   :  return error : 0 - OK,
	                   1 - opening problem,
                           2 - file longer then array.
        fname : file name.
        size  : on input - max size of the arr or 0 for any length,
                on output- real size of the file (and arr in bytes).
        arr   : returned file as array.

  History:  A.Jesmanowicz, MCW 1991
*/
void read_iqm( long int *ier, char fname[], long int *size, char arr[] )
{
        int     isize = *size;
        int     fp;                             /* file descriptor  */
        struct stat file_stat;                  /* status structure */
 
        {       int i=0;                        /* cut junk from the name */
                while(fname[i] != 32 && fname[i]) i++;
                fname[i] = 0;
        }
 
        if ((fp = open(fname, O_RDONLY)) <= 0)  /* file must exist */
        {       *ier=1;                         /* or ier = 1.     */
		printf("File %s does not exist!\n",fname);
		return;
        }
 
        fstat(fp, &file_stat);                  /* get file size in bytes   */
 
        if(file_stat.st_size > isize && isize)  /* file can not be too long */
        {       *ier=2;                         /* or ier = 2.              */
		printf("File %s is too long!\n",fname);
		exit(0);
        }
        *size =  file_stat.st_size;             /* return file size */
 
        read(fp, arr, file_stat.st_size);       /* read whole file  */
        close(fp);
        *ier=0;                                 /* no error : ier=0 */
        return;
}
/**<write_iqm.c>***********************************************
 
    Title: WRITE I/Q channels
 
  Purpose: Write complex data

     Call: void write_iqm(ier,fname,size,arr)
             long int *ier;
	     long int *size;
	     char fname[],arr[];

   ier   :  return error : 0 - OK,
	                   1 - opening problem,
                           2 - file longer then array.
        fname : file name.
        size  : on input - max size of the arr or 0 for any length,
                on output- real size of the file (and arr in bytes).
        arr   : returned file as array.

  History:  A.Jesmanowicz, MCW 1991
            lrf 24jun96  fixed overwrite bug
*/
void write_iqm( long int *ier, char fname[], long int *size, char arr[] )
{
	int	isize = *size;
	int	fp;				/* file descriptor  */

	{	int i=0;			/* cut junk from the name */
		while(fname[i] != 32 && fname[i]) i++;
		fname[i] = 0;
	}

	if(isize < 0)				/* size has to be real */
	{	*ier=2;				/* or ier = 2.	       */
		return;
	}

	if ((fp = open(fname, O_WRONLY|O_CREAT|O_TRUNC,0777)) <= 0)
	{	*ier=1;
		printf("Error in open = %i\n",(int)*ier);
		return;				/* or ier = 1.	   */
	}

	*ier=0;					/* no error : ier=0 */
	write(fp, arr, isize);  		/* write whole file */
	close(fp);
	return;  
}
