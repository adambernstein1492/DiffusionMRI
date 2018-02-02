/* excite.h */

double *minfinity_rot(double *M,double theta,double phi,double xi,double we,double t1,double t2) ;
double *propagate_rot( double *M,double **R,double **E,
		       double w1,double p1,double dw,double t1,double t2,double dt) ;
double **rotation_matrix_fm( double **R, double phi ) ;
double **rotation_matrix_rot_3d(double **R,double theta,double phi,double xi) ;
double **relaxation_matrix_rot_3d( double **E, double theta, double phi, 
				   double xi, double we, double t1, double t2) ;


