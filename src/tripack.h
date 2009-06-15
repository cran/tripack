
#include <R.h>

int F77_NAME(addcst)(int *ncc, int *lcc, int *n, 
		     double *x, double *y, int *lwk, int *iwk, 
		     int *list, int *lptr, int *lend, int *ier);

int F77_NAME(troutq)(int *ncc, int *lcc, int *n, 
		     double *x, double *y, int *list, int *lptr, 
		     int *lend, int *lout, int *nnabs, int *nptr, 
		     int *nptr1, int *nabor, int *nbnos, int *na, 
		     int *nb, int *nt);

int F77_NAME(bnodes)(int *n, int *list, int *lptr, int *lend, 
		     int *nodes, int *nb, int *na, int *nt); 

int F77_NAME(inhull)(double *xp, double *yp, int *np, 
		     double *x, double *y, int *n, int *list, 
		     int *lptr, int *lend, int *inh);

int F77_NAME(onhull)(double *xp, double *yp, int *np, 
		     double *x, double *y, int *n, int *list, 
		     int *lptr, int *lend, int *onh, double *eps);

int F77_NAME(trlist)(int *ncc, int *lcc, int *n, int *list, 
		     int *lptr, int *lend, int *nrow, int *nt,
		     int *ltri, int *lct, int *ier);

int F77_NAME(trfind)(int *nst, double *px, double *py, int *n, 
		     double *x, double *y, int *list, int *lptr, 
		     int *lend, int *i1, int *i2, int *i3);

int F77_NAME(trmesh)(int *n, double *x, double *y, int *list, 
		     int *lptr, int *lend, int *lnew, int *near, 
		     int *next, double *dist, int *ier); 

int F77_NAME(voronoi)(int *ncc, int *lcc, int *n, 
		      double *x, double *y, int *list, int *lptr, 
		      int *lend, int *nt, double *lccc, int *iccc, 
		      int *lct, int *ltri, int *ier);
