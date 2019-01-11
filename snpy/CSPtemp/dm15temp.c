/* dm15temp.c:  driver program for GLoEs (Gaussial Local Estimation) to 
 * smooth the 2-d surface of time-dm15-flux  lightcurve data.  
 *
 * Author:  Chris Burns (cburns@ociw.edu)
 *
 * The problem:  you've got a set of low-z lightcurves in certain filters.  You
 * have analyzed each lightcurve to determine the time of maximum, the maximum
 * magnitude, and the decline rate (dm15).  Now, each observation can be placed
 * on a grid (x -> epoch, y-> dm15) and the flux (magnitude) defines a surface
 * on this grid.  We now wish to sample this surface (defined by heterogenously
 * spaced points) at any arbitrary (epoch,dm15) point.  Typically, this is so 
 * that you can fit another lightcurve (high-z or without pre-max data) 
 * using least squares.
 *
 * The solution:  Use GLoEs (Gaussian Local Estimation) to smooth the 2D
 * surface and interpolate to any given point you wish.  It does this by fitting
 * a 2nd degree polynomial centered at the point of interest.  The data points
 * are weighted by both their internal errors as well as a sliding 2D
 * elliptical Gaussian of variable width (hence the "Gaussian" in the name).  In
 * this way, only points close to the interpolating point are used (hence the
 * "local" in the name).  See gloes.c for more details about the algorithm.
 *
 * This program expects to find a file named "templates.dat" in which you have
 * listed the locations of the template lightcurve data.  The format of
 * templates.dat should be:
 * SN-name ; filter ; redshift ; dm15 ; T(Bmax) ; Mmax ; filename
 * (ie, semi-colon separated).  SN-name really isn't used, it's just for your
 * own sanity.  Here's an example:
 *
 * SN2004gf ; B ; 0.01 ; 1.1 ; 452.5 ; 16.5 ; SN2004gf_B.dat
 * SN2004gf ; V ; 0.01 ; 1.1 ; 452.5 ; 15.9 ; SN2004gf_V.dat
 * SN2004gf ; g ; 0.01 ; 1.1 ; 452.5 ; 15.8 ; SN2004gf_g_prime.dat
 * SN2004gf ; r ; 0.01 ; 1.1 ; 452.5 ; 15.4 ; SN2004gf_r_prime.dat
 * ....
 *
 * Then, each datafile listed in this file has the following format:
 * time magnitude magnitude-error
 *
 * The time should be in days.  Tmax listed in templates.dat will be subtracted
 * to get the epoch.  The Mmax from templates.dat will be subtracted from the
 * magnitudes so that Mmax = 0 for each template (or M(Tmax) = 0 if you choose
 * Mmax accordingly).
 * 
 * This program will read in templates.dat and then the files listed therein,
 * filling in the arrays needed to define the surface.  It reads the value of
 * dm15 from the command line and then samples the surface on the interval 
 * [-10,70] in one-day increments.  The interpolated values and errors are
 * printed to STDOUT. 
 *
 * The width of the window function is a function of position on the grid.
 * Currently, a fixed formula is used (see sigmax() and sigmay()), since
 * automatic algorithms in 2D are non-trivial.  These formulas are purely
 * trial and error:  they seem to work well.  As data are added, they will
 * likely need to be updated to improve efficiency (or experiment for
 * yourself).  But generically, one can say that there is very little
 * curvature in the dm15-direction (hence I use a fixed width), and high
 * curvature near max in epoch-direction, which decreases away from max,
 * hence I use a linearly increasing window.
 *
 * You need the GSL libraries and include files.  Compile using:
 *
 * gcc -o dm15temp dm15temp.c gloes.c -lm -lgsl -lgslcblas -I{Idir} -L{Ldir}
 *
 * where {Idir} and {Ldir} should be replaced with the location of the gsl
 * include files and libraries, respectfully, if they are not in the standard
 * search path.  Then, simply call the program like:
 *
 * dm15temp dm15 dm15 dm15 > data
 *
 * You can list any number of values for dm15, they will be output
 * sequencially to STDOUT.  The format is described in the output header.  Note that
 * the lightcurves are in flux units and are normalized to 1.0 at the maximum of
 * the _individual lightcurves_.  So, at t(Bmax), the SN will not, in general, have
 * 0 color.  It will, however, have all pseudo-colors at max (e.g., Bmax-Vmax) = 0
 * so it's up to you to impose any colors yourself (e.g., Phillips et al, 1999).
 *
 * */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_min.h>
#include "gloes.h"

#define n_filters 16
#define max_points 5000

 /* global arrays to contain the data */
double x[n_filters][max_points], y[n_filters][max_points], 
       z[n_filters][max_points], wz[n_filters][max_points];
int np[n_filters];
char name[n_filters][max_points][80];

double sigx0;
double sigy0;
double xscale;
double sigxmax;

/* function for returning the width of the window function in epoch as a function
 * of epoch and dm15.  Right now, just a linearly increasing window with a cap */
double  sigmax(double x, double y) {
   double sigma = sigx0 + xscale*abs(x);
   if (sigma > sigxmax) sigma = sigxmax;
   return(sigma);
}

/* function for returning the width of the window function in dm15 as a function
 * of epoch and dm15.  Right now, just a constant value */
double sigmay(double x, double y) {
   double sigma = sigy0;
   return(sigma);
}

int load_data(char *path) {
   FILE *fp1, *fp2;
   int i, j, nchars;
   char filename[1024];
   char filename2[256];
   char line[256];
   char SN[80];
   char filter[15];
   double dm15, Tmax, Mmax, zed;
   double epoch, mag, emag;

   nchars = snprintf(filename, 1024, "%s/%s", path, "templates.dat");
   if (nchars > 1024) {
      fprintf(stderr, "Error:  path name for templates.dat is too long\n");
      return(-1);
   }

   if (! (fp1 = fopen(filename, "r"))) { 
      fprintf(stderr, "%s not found!\n", filename);
      exit(1);
   }
   /* Initialize arrays */
   for (i = 0 ; i < n_filters ; ++i) {
      np[i] = 0;
      for(j = 0 ; j < max_points ; ++j) {
         x[i][j] = 0.0;
         y[i][j] = 0.0;
         z[i][j] = 0.0;
         wz[i][j] = 0.0;
      }
   }

   /* Read in template information */
   while(fgets(line, 256, fp1) != NULL) {
      if(line[0] == '#') continue;
      line[255] = '\0';
      strcpy(SN, strtok(line, ";"));
      sscanf(strtok(NULL, ";"), "%s", filter);
      zed = atof(strtok(NULL, ";"));
      dm15 = atof(strtok(NULL, ";"));
      Tmax = atof(strtok(NULL, ";"));
      Mmax = atof(strtok(NULL, ";"));
      sscanf(strtok(NULL, ";"), "%s", filename2);
      nchars = snprintf(filename, 1024, "%s/%s", path, filename2);
      if (nchars > 1024) {
         fprintf(stderr, "Error:  path name for template is too long\n");
         return(-1);
      }
      if ( ! ( fp2 = fopen(filename, "r") ) ) {
         fprintf(stderr, "Could not find file '%s'\n", filename);
         exit(1);
      }
      if (strcmp(filter, "B") == 0) {
         i = 0;
      } else if (strcmp(filter, "V") == 0) {
         i = 1;
      } else if (strcmp(filter, "u") == 0) {
         i = 2;
      } else if (strcmp(filter, "g") == 0) {
         i = 3;
      } else if (strcmp(filter, "r") == 0) {
         i = 4;
      } else if (strcmp(filter, "i") == 0) {
         i = 5;
      } else if (strcmp(filter, "Y") == 0) {
         i = 6;
      } else if (strcmp(filter, "J") == 0) {
         i = 7;
      } else if (strcmp(filter, "H") == 0) {
         i = 8;
      } else if (strcmp(filter, "K") == 0) {
         i = 9;
      } else {
         fprintf(stderr, "Error, filter %s not recognized\n", filter);
         exit(1);
      }

      while( fscanf(fp2, "%lf %lf %lf", &epoch, &mag, &emag) != EOF) {
         x[i][np[i]] = (epoch - Tmax)/(1.0+zed);
         y[i][np[i]] = dm15;
         z[i][np[i]] = pow(10, -0.4*(mag - Mmax));
         wz[i][np[i]] = emag > 0 ? 1.0857/(z[i][np[i]]*emag) : 0.0;  /* check for 1/0 */
         nchars = snprintf(name[i][np[i]], 80, "%s", SN);
         np[i] += 1;
      }
      fclose(fp2);
   }
   fclose(fp1);
   return(0);
}

struct my_params { double dm15; int filter; } ;

double template(double t, void *p) {
   struct my_params *params = (struct my_params * ) p;
   double sigx, sigy;
   double dm15 = (params->dm15);
   int filter = (params->filter);
   int res;
   double mag, emag;
   sigx = sigmax(t, dm15);
   sigy = sigmay(t, dm15);
   res = gloes2D_n_sigma(x[filter], y[filter], z[filter], wz[filter], np[filter],
         &t, &dm15, 1, &sigx, &sigy, &mag, &emag);
   return (-1.0*mag);
}


int dm15temp(int filter, double dm15, double *t, int size, double *mag,
      double *emag, int normalize, double arg1, double arg2, double arg3,
      double arg4) {
   int i;
   int res;
   double *sigx, *sigy;
   struct my_params p;
   double *dm15s;
   const gsl_min_fminimizer_type *T;
   gsl_min_fminimizer *s;
   int iter=0, maxiter=100;
   int status;
   double m, a, b, max;
   gsl_function func;

   sigx0 = arg1;
   sigy0 = arg4;
   xscale = arg2;
   sigxmax = arg3;

   sigx = (double *) malloc(size*sizeof(double));
   sigy = (double *) malloc(size*sizeof(double));
   dm15s = (double *) malloc(size*sizeof(double));
   for (i = 0 ; i < size ; ++i ) {
      sigx[i] = sigmax(t[i], dm15);
      sigy[i] = sigmay(t[i], dm15);
      dm15s[i] = dm15;
   }

   res = gloes2D_n_sigma(x[filter], y[filter], z[filter], wz[filter], 
         np[filter], t, dm15s, size, sigx, sigy, mag, emag);
   if (normalize) {
      gsl_set_error_handler_off();
      /* Normalize the lightcurves ...  find the maximum */
      func.function = &template;
      func.params = &p;
      T = gsl_min_fminimizer_brent;
      s = gsl_min_fminimizer_alloc(T);
      p.dm15 = dm15;
      p.filter = filter;
      status = gsl_min_fminimizer_set(s, &func, 0.0, -10.0, 20.0);
      if (status == GSL_EINVAL) {
         /* fprintf(stderr, "Normalization Warning:  gsl_min_fminimizer failed.\n"); */
         gsl_min_fminimizer_free(s);
         free(sigx);
         free(sigy);
         free(dm15s);
         return(5);
      }

      do {
         iter++;
         status = gsl_min_fminimizer_iterate(s);
         m = gsl_min_fminimizer_x_minimum(s);
         max = -1.0*gsl_min_fminimizer_f_minimum(s);
         a = gsl_min_fminimizer_x_lower(s);
         b = gsl_min_fminimizer_x_upper(s);
         status = gsl_min_test_interval(a, b, 0.0001, 0.0);
      } while (status == GSL_CONTINUE && iter < maxiter);
      for (i = 0 ; i < size ; ++i ) {
         mag[i] = mag[i]/max;
         emag[i] = emag[i]/max;
      }
      gsl_min_fminimizer_free(s);
   }
   free(sigx);
   free(sigy);
   free(dm15s);
   return(res);
}

    
int main(int argc, char **argv) {
   double xp[100], zp[n_filters][100], ezp[n_filters][100];
   int i,j,res;
   double dm15;

   if (argc < 2) {
      fprintf(stderr, "Usage:  dm15temp dm15 dm15 dm15 ...\n");
      fprintf(stderr, "       where dm15 are the (possible many) decline rates\n");
      exit(1);
   }
   load_data(".");

   for (j = 1 ; j < argc ; ++j ) {
      /* setup the evaluation points, result arrays and window widths */
      dm15 = atof(argv[j]);
      for (i = 0 ; i < 81 ; ++i) {
         xp[i] = -10.0 + i;
      }
 
      for (i = 0 ; i < 9 ; ++i) {
         res = dm15temp(i, dm15, xp, 81, zp[i], ezp[i], 1, 3.0, 0.1, 10.0, 0.3);
      }

      printf("# Lightcurve generated for dm15 = %f\n", dm15);
      printf("# col(1) = epoch (days)\n");
      printf("# col(2) = B mag (mag)\n");
      printf("# col(3) = error in B mag (mag)\n");
      printf("# col(4) = V mag (mag)\n");
      printf("# col(5) = error in V mag (mag)\n");
      printf("# col(6) = u mag (mag)\n");
      printf("# col(7) = error in u mag (mag)\n");
      printf("# col(8) = g mag (mag)\n");
      printf("# col(9) = error in g mag (mag)\n");
      printf("# col(10) = r mag (mag)\n");
      printf("# col(11) = error in r mag (mag)\n");
      printf("# col(12) = i mag (mag)\n");
      printf("# col(13) = error in i mag (mag)\n");
      printf("# col(12) = Y mag (mag)\n");
      printf("# col(13) = error in Y mag (mag)\n");
      printf("# col(12) = J mag (mag)\n");
      printf("# col(13) = error in J mag (mag)\n");
      printf("# col(12) = H mag (mag)\n");
      printf("# col(13) = error in H mag (mag)\n");
      for (i = 0 ; i < 81 ; ++i) {
         printf("%4.1f  %7.4f %6.4f %7.4f %6.4f %7.4f %6.4f %7.4f %6.4f %7.4f %6.4f %7.4f %6.4f %7.4f %6.4f %7.4f %6.4f %7.4f %6.4f\n", 
               xp[i], zp[0][i], ezp[0][i], zp[1][i], ezp[1][i],zp[2][i], ezp[2][i],zp[3][i], ezp[3][i],zp[4][i],
               ezp[4][i],zp[5][i], ezp[5][i], zp[6][i], ezp[6][i],zp[7][i], ezp[7][i],zp[8][i], ezp[8][i]);
      }
      printf("\n");
   }
   return(0);
}
