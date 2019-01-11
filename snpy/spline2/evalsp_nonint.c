/*
 *      Evalsp
 *
 *      Program for evaluating, analyzing, and manipulating a spline of which
 *   the so-called beta-spline representation (BS representation) is
 *      given. The BS-representation consists of the order k, the number 
 *   of intervals l, the values of the l+1 breakpoints, and the l+k-1 
 *   beta-spline coefficients.
 *   
 *   This program is intended to be used after running "spline2", which
 *   determines such a spline as a least-squares approximation to
 *   data points, and writes the above BS-representation to a file 
 *   named "splrep".
 *
 *  16 Jan 2002: Default number of output points increased from 250 to 2000.
 *   5 Feb 2002: Cosmetic changes to screen and evlsum layout.
 *
 ****   Syntax:
 *
 *   evalsp [filename1] [-x xbegin xend xstep] [-f xfile] [-v value]
 *          [-e] [-i] [-a] [-F sbegin send sstep] [-h] [-b value]
 *
 ****   Command line options:
 *
 *   [filename1]      Name of file containing BS-representation,
 *            usually called "splrep".
 *            (Default=stdin).
 *
 *   [-x xbegin xend xstep]   Values in which the spline and its 
 *            derivatives are to be evaluated.
 *            This is essentially an interpolation operation.
 *            (Default: in 2000 equidistant steps between
 *            the first and last breakpoints, inclusive.
 *            This is convenient for plotting).
 *
 *            Note: the spline and its derivatives evaluated
 *            in the x-values of the original datapoints
 *            can be found in file "splres". This file is
 *            created when "spline2" is used with the -F
 *            option.
 *
 *   [-f xfile]      Name of file containing x-values (in
 *            increasing order) in which the spline and its
 *            derivatives are to be evaluated. This option
 *            is an alternative to the -x option. It allows
 *            the spline to be evaluated in x-values that
 *            belong, for instance, to another dataset.
 *
 *   [-v value]      Generates the x-value(s) in which the spline
 *            is equal to "value" (inverse interpolation,
 *            using Brent's method).
 *
 *   [-e]         Generates the x- and y-values of the minima
 *            and maxima of the spline.
 *
 *   [-i]              Generates the x-, y- and dy/dx-values of the
 *            inflection points of the spline.
 *
 *   [-a]         Calculates the integral from the first to
 *            the last breakpoint.
 *
 *   [-F sbegin send sstep]  Evaluates the Fourier transform of the
 *            spline from sbegin to send in steps of sstep.
 *            Note: a direct fft is often to be preferred.
 *
 *   [-h]         Generates the x-values in which the spline
 *            reaches the "half-height" value, which
 *            is assumed to be the value halfway the
 *            highest maximum and the avarage value
 *            of all other minima and maxima. This is
 *            a very primitive method to determine the width
 *            of a single peak (on a possibly fluctuating
 *            background), but it works sometimes.
 *            Implies -e (and later: -v).
 *
 *   [-b value]      Specifies a baseline (background) value
 *            instead of having it calculated by the -h
 *            option.
 *
 *   Note concerning options -v, -e and -i: if the range between two 
 *   subsequent x-values contains more than one root, no root is found 
 *   in the case of an even number, and only one root is found in the case
 *   of an odd number of roots.
 *
 ****   Output:
 *
 *      File "evlres", containing x(i), S(x(i)), and all derivatives
 *      dj(S(i))/dxj; from j=1 to j=k-1. One line for each
 *      x-value.
 *
 *   File "evlsum", containing the results of options -v, -e, -i, -h and -a.
 *      This information is also written to stdout (the terminal).
 *
 *   File "evlfour", containing the values of the Fourier transform 
 *      as follows: s, Re(Ft), Im(Ft). One line for each
 *              s-value.
 *
 ****   Info:
 *
 *   Authors and owners of this program are:
 *      Barend J. Thijsse and Mark A. Hollanders,
 *              Computational Materials Science group (Com,ma,s)
 *      Physical Materials Science Division,
 *      Laboratory of Materials Science, 
 *      Delft University of Technology,
 *      Rotterdamseweg 137, 2628 AL  Delft,
 *         Phone: +31 15 278 2221
 *         Fax:   +31 15 278 6730
 *         E-mail: thijsse@stm.tudelft.nl
 *         URL: http://dutsm183.stm.tudelft.nl
 *
 */

/* ---INCLUDES--- */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

/* ---DEFINES (NUMBERS)--- */
#define NMAX 10001
#define LMAX 1001
#define KMAX 24
#define ROOTMAX 100

/* ---ANNOUNCE FUNCTIONS--- */ 
int    l2knts(),
       bsplpp(),
       bsplvb(),
       interv();
double ppvalu(),
       ppigr(),
       zbrent();

/* ---GLOBAL VARIABLES (are zero-initialized)--- */
int    k, km1, l, n, left, jbsp, ihi, ilo = 1, hhon, nextremes, nomaxmax=1;
int    firsttime=1, very_uncertain;
double brek[LMAX], bcoef[LMAX+KMAX], coef[KMAX][LMAX], ppvalue,
       deltal[KMAX], deltar[KMAX]; 
double halfheight, extrsum, maxmax, zerolevel;
FILE   *fpw, *fps, *fpf, *fpr, *fpx;



/* ---MAIN--- */
int eval_extrema(int in_k, int in_l, double *in_t, double *in_c, double **xextr, double **yextr, int **signs,
      int *n_extr){
   register int i, j;
   int          xflag, ntau, m, fileind, val, ext, inf, are, fou;
   int        background;
   double       xbegin, xend, xstep, y, x, t[LMAX+KMAX], b, bb, a[KMAX],
           scrtch[KMAX][KMAX], c, a0, a1, a2, yder, ymax, step,
           xv[NMAX], value, smin, smax, sstep, pi, s, fourr, fouri,
           tpis, pref, argb, arga, bg;
   char        str[100];
   double xmaxs[ROOTMAX], ymaxs[ROOTMAX];
   int maxflags[ROOTMAX];

   /* ---INITIALIZATION---*/
     xflag = fileind = val = ext = inf = are = fou = 0;
     background = 0;

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   if (k <= 2)
   {
      fprintf(stderr, " spline order too low for finding extrema\n");
      return(1);
   }
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   ntau = 2001;
   xbegin = brek[1];
   xend = brek[l+1];
   xstep = (xend - xbegin)/2000;

   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);
   nextremes = 0;
   for (i = 1; i <= ntau; i++)
   {
      b = xbegin + (i-1)*xstep;
      a[0] = ppvalu(brek,coef,&l,&k,&b,0);
      a[1] = ppvalu(brek,coef,&l,&k,&b,1);
      a[2] = ppvalu(brek,coef,&l,&k,&b,2);

      if (i > 1 && a[1]*a1 <= 0. && a1 != 0 && a[2] != 0)
      {
         if (a[1] == 0.) x = b;
         else x = zbrent(bb,b,1,0);
         ymax = ppvalu(brek,coef,&l,&k,&x,0);
         if (a1 > 0.)
         {
         xmaxs[nextremes] = x;
         ymaxs[nextremes] = ymax;
         maxflags[nextremes] = -1;
         nextremes++;
         }
         else
         {
         xmaxs[nextremes] = x;
         ymaxs[nextremes] = ymax;
         maxflags[nextremes] = 1;
         }
      }
      a0 = a[0];
      a1 = a[1];
      a2 = a[2];
      bb = b;
   }
   if (! (*xextr = (double *)malloc(*n_extr*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_y in eval_inflect\n");
      return(1);
   }
   if (! (*yextr = (double *)malloc(*n_extr*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_y in eval_inflect\n");
      return(1);
   }
   if (! (*signs = (int *)malloc(*n_extr*sizeof(int)))) {
      fprintf(stderr, "Bad MALLOC for out_y in eval_inflect\n");
      return(1);
   }
   *n_extr = nextremes;
   for (i = 0 ; i < nextremes ; ++i) {
      *xextr[i] = xmaxs[i];
      *yextr[i] = ymaxs[i];
      *signs[i] = maxflags[i];
   }
   return(0);
}



int eval_inflect(int in_k, int in_l, double *in_t, double *in_c, double **out_y, double **out_dy, 
      int *n_inflect){
   register int i, j;
   int          xflag, ntau, m, fileind, val, ext, inf, are, fou;
   int        background;
   double       xbegin, xend, xstep, y, x, t[LMAX+KMAX], b, bb, a[KMAX],
           scrtch[KMAX][KMAX], c, a0, a1, a2, yder, ymax, step,
           xv[NMAX], value, smin, smax, sstep, pi, s, fourr, fouri,
           tpis, pref, argb, arga, bg;
   double ymaxs[ROOTMAX], yders[ROOTMAX];
   char        str[100];

   /* ---INITIALIZATION---*/
     xflag = fileind = val = ext = inf = are = fou = 0;
     background = 0;

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   if (k <= 3)
   {
      fprintf(stderr, " spline order too low for inflection evaluation\n");
      return(1);
   }
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   ntau = 2001;
   xbegin = brek[1];
   xend = brek[l+1];
   xstep = (xend - xbegin)/2000.0;

   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);
   *n_inflect = 0;
   for (i = 1; i <= ntau; i++)
   {
      b = xbegin + (i-1)*xstep;
      a[0] = ppvalu(brek,coef,&l,&k,&b,0);
      a[1] = ppvalu(brek,coef,&l,&k,&b,1);
      a[2] = ppvalu(brek,coef,&l,&k,&b,2);
      a[3] = ppvalu(brek,coef,&l,&k,&b,3);

      if (i > 1 && a[2]*a2 <= 0. && a2 != 0 && a[3] != 0)
      {
         if (a[2] == 0.) x = b;
         else x = zbrent(bb,b,2,0);
         ymax = ppvalu(brek,coef,&l,&k,&x,0);
         yder = ppvalu(brek,coef,&l,&k,&x,1);
         if(yder>0.0){
            ymaxs[*n_inflect] = ymax;
            yders[*n_inflect] = yder;
            *n_inflect += 1;
         }
         if(yder<0.0){
            ymaxs[*n_inflect] = ymax;
            yders[*n_inflect] = yder;
            *n_inflect += 1;
         }
         if(yder==0.0){
            ymaxs[*n_inflect] = ymax;
            yders[*n_inflect] = yder;
            *n_inflect += 1;
         }
      }
      a0 = a[0];
      a1 = a[1];
      a2 = a[2];
      bb = b;
   }
   if (! (*out_y = (double *)malloc(*n_inflect*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_y in eval_inflect\n");
      return(1);
   }
   if (! (*out_dy = (double *)malloc(*n_inflect*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_dy in eval_inflect\n");
      return(1);
   }
   for (i = 0 ; i < *n_inflect ; ++i ) {
      *out_y[i] = ymaxs[i];
      *out_dy[i] = yders[i];
   }
   return(0);
}


double eval_integ(double x0, double x1, int in_k, int in_l, double *in_t, double *in_c) {
   register int i, j;
   int          xflag, ntau, m, fileind, val, ext, inf, are, fou;
   int        background;
   double       xbegin, xend, xstep, y, x, t[LMAX+KMAX], b, bb, a[KMAX],
           scrtch[KMAX][KMAX], c, a0, a1, a2, yder, ymax, step,
           xv[NMAX], value, smin, smax, sstep, pi, s, fourr, fouri,
           tpis, pref, argb, arga, bg;
   char        str[100];

   /* ---INITIALIZATION---*/
     xflag = fileind = val = ext = inf = are = fou = 0;
     background = 0;

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   xbegin = brek[1];
   xend = brek[l+1];
   if (x0 < xbegin || x1 > xend) {
      fprintf(stderr, "Error:  integration limits beyond spline definition\n");
      return(1);
   }
   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);

   y = ppigr(brek,coef,&l,&k,&xbegin,&xend);
   return(y);
}


int eval_x(double value, int in_k, int in_l, double *in_t, double *in_c,double **ret_roots, int *n_roots) {
   register int i, j;
   int          ntau;
   double  xbegin, xend, xstep, y, x, t[LMAX+KMAX], b, bb, a[KMAX],
           scrtch[KMAX][KMAX], a0, yder, ymax, step,
           xv[NMAX], pi, s,
           tpis, pref;
   double roots[ROOTMAX];

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   if (k <= 1)
   {
      fprintf(stderr, "spline order too low for reverse evaluation\n");
      return(1);
   }
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   ntau = 2001;
   xbegin = brek[1];
   xend = brek[l+1];
   xstep = (xend - xbegin)/2000;

   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);
   *n_roots = 0;
   for (i = 1; i <= ntau; i++)
   {
      b = xbegin + (i-1)*xstep;
      a[0] = ppvalu(brek,coef,&l,&k,&b,0);

      if (i > 1 && (a[0]-value)*(a0-value) <= 0. && a0 != 0)
      {
         if (a[0]-value == 0.) x = b;
         else x = zbrent(bb,b,0,value);
         roots[*n_roots] = x;
         *n_roots += 1;
      }
      a0 = a[0];
      bb = b;
   }
   if (! (*ret_roots = (double *)malloc(*n_roots*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for ret_roots in eval_x\n");
      return(1);
   }
   for(i = 0 ; i < *n_roots ; ++i) {
      *ret_roots[i] = roots[i];
   }
   return(0);
}


int evalsp(double *in_x, int num_points, int in_k, int in_l, double *in_t, double *in_c, int deriv, 
      double **out_y) {
   /* Evaluate the spline.  *x is the input num_points x-values.  If values=1, then treat x-values as
    * y-values and do a reverse-interpolation.  If deriv>0, return the deriv-th derivative of the
    * interpolation (up to k-1) */
   /* ---DECLARATIONS (note: auto-arrays can't be initialized)--- */

   register int i, j;
   int          xflag, ntau, m, fileind, val, ext, inf, are, fou;
   int        background;
   double       xbegin, xend, xstep, y, x, t[LMAX+KMAX], b, bb, a[KMAX],
           scrtch[KMAX][KMAX], c, a0, a1, a2, yder, ymax, step,
           xv[NMAX], value, smin, smax, sstep, pi, s, fourr, fouri,
           tpis, pref, argb, arga, bg;
   char        str[100];

   /* ---INITIALIZATION---*/
     xflag = fileind = val = ext = inf = are = fou = 0;
     background = 0;

   /* ---PROCEED---*/

   k = in_k;
   l = in_l;
   for (i = 1; i <= l+1; i++) brek[i] = in_t[i-1];
   km1 = k-1;
   for (i = 1; i <= l+km1; i++) bcoef[i] = in_c[i-1];
   for (i = 1 ; i <= num_points ; ++i) xv[i] = in_x[i-1];
   ntau = num_points;
   xbegin = xv[1];
   xend = xv[ntau];

   l2knts(brek,&l,&k,t,&n);
   bsplpp(t,bcoef,&n,&k,scrtch,brek,coef,&l);

   if( ! (*out_y = (double *) malloc(num_points*sizeof(double)))) {
      fprintf(stderr, "Bad MALLOC for out_y in evalsp\n");
      return(1);
   }

   for (i = 1; i <= ntau; i++)
   {
      b = xv[i];
      *out_y[i-1] = ppvalu(brek,coef,&l,&k,&b,deriv);
   }

   return(0);
   } /* --end main-- */

/*------------------------------------------------------------------*/
double ppigr(ara,dara,iptr,jptr,xptr,yptr)
   int     *iptr, *jptr;
   double  *xptr, *yptr, ara[], dara[][LMAX];

   /*  calculates integral from *xptr to *yptr  */
{
   int     i, j, right, ndummy;
   double  h, aa, bb, ppintgr;
   ppintgr = 0.;
   interv(ara,iptr,xptr,&left,&ndummy);
   h = *xptr-ara[left];
   bb = 0.;
   for (j = *jptr; j >= 1; j--) bb = (bb+dara[j][left])*h/j;
   interv(ara,iptr,yptr,&right,&ndummy);
   for (i = left; i < right; i++)
   {
      aa = 0.;
      h = ara[i+1]-ara[i];
      for (j = *jptr; j >= 1; j--)
      {
         aa = (aa+dara[j][i])*h/j;
      }
      ppintgr += aa;
   }
   h = *yptr-ara[right];
   aa = 0.;
   for (j = *jptr; j >= 1; j--) aa = (aa+dara[j][right])*h/j;
   ppintgr = ppintgr+aa-bb;
   return(ppintgr);
}

/*------------------------------------------------------------------*/
int l2knts(ara,iptr,jptr,arb,kptr)
   int   *iptr, *jptr, *kptr;
   double  *ara, *arb;  

   /*  breakpoints to knots  */
{
   int   i;
   for (i = 1; i <= km1; i++) arb[i] = ara[1];
   for (i = 1; i <= *iptr; i++) arb[km1+i] = ara[i];
   n = km1+(*iptr);
   for (i = 1; i <= *jptr; i++) arb[*kptr+i] = ara[*iptr+1];
}

/*------------------------------------------------------------------*/
int bsplpp(ara,arb,iptr,jptr,dara,arc,darb,kptr)
   int   *iptr, *jptr, *kptr;
   double  ara[], arb[], arc[], dara[][KMAX], darb[][LMAX];

   /*  converts spline to piecewise polynomial representation  */
{
   int   lsofar, j, i, jp1, kmj;
   double  diff, sum, biatx[KMAX];
   arc[1] = ara[*jptr];
   lsofar = 0;
   for (left = *jptr; left <= *iptr; left++)
   {
      if(ara[left+1] != ara[left])
      {
         lsofar++;
         arc[lsofar+1] = ara[left+1];
         if (*jptr <= 1) darb[1][lsofar] = arb[left];
         else
         {
            for (i = 1; i <= *jptr; i++)
            {
               dara[i][1] = arb[left-*jptr+i];
            }
            for (jp1 = 2; jp1 <= *jptr; jp1++)
            {
               j = jp1-1;
               kmj = k-j;
               for(i = 1; i <= kmj; i++)
               {
               diff = ara[left+i]-ara[left+i-kmj];
                  if (diff > 0.)
                  {
         dara[i][jp1] = ((dara[i+1][j]-dara[i][j])/diff)*kmj;
                  }
               }
            }
            bsplvb(ara,1,1,&ara[left],&left,biatx);
            darb[*jptr][lsofar] = dara[1][*jptr];
            for(jp1 = 2; jp1 <= *jptr; jp1++)
            {
            bsplvb(ara,jp1,2,&ara[left],&left,biatx);
               kmj = k+1-jp1;
               sum = 0.;
               for(i = 1; i <=jp1; i++)
               {
                  sum += biatx[i]*dara[i][kmj];
                  darb[kmj][lsofar] = sum;
               }
            }
         }
      }
   }
   *kptr = lsofar;
}

/*------------------------------------------------------------------*/
int bsplvb(ara,jhigh,index,xptr,iptr,arb)
   int   jhigh, index, *iptr;
   double  ara[], arb[], *xptr;

   /*  calculates all nonzero beta-splines at *xptr  */
{
   int   jp1, i;
   double  saved, term;
   if (index == 1)
   {
      jbsp = 1;
      arb[1] = 1.;
   }
   while (jbsp < jhigh)
   {
      jp1 = jbsp+1;
      deltar[jbsp] = ara[*iptr+jbsp]-*xptr;
      deltal[jbsp] = (*xptr)-ara[*iptr+1-jbsp];
      saved = 0.;
      for (i = 1; i <= jbsp; i++)
      {
         term = arb[i]/(deltar[i]+deltal[jp1-i]);
         arb[i] = saved+deltar[i]*term;
         saved = deltal[jp1-i]*term;
      }
      arb[jp1] = saved;
      jbsp++;
   }
}

/*------------------------------------------------------------------*/
double ppvalu(ara,dara,iptr,jptr,xptr,jderiv)
   int   *iptr, *jptr, jderiv;
   double  *xptr, ara[], dara[][LMAX];

   /*  evaluates the jderiv-th derivative of a pp-function  */
{
   int   fmmjdr, i, ndummy, m;
   double  h;
   ppvalue = 0.;
   fmmjdr = *jptr-jderiv;
   if (fmmjdr > 0)
   {
      interv(ara,iptr,xptr,&i,&ndummy);
      h = *xptr-ara[i];
      for (m = *jptr; m >= jderiv+1; m--)
      {
         ppvalue = (ppvalue/fmmjdr)*h+dara[m][i];
         fmmjdr--;
      }
   }
   return(ppvalue);
}

/*------------------------------------------------------------------*/
int interv(ara,iptr,xptr,jptr,kptr)
   int   *iptr, *jptr, *kptr;
   double  *xptr, ara[];

   /*  locates a point within an increasing sequence of points  */
{
   int   istep, middle, ilos;
   *kptr = 10;
   ihi = ilo+1;
   if (ihi >= *iptr)
   {
      if (*xptr >= ara[*iptr]) *kptr = 1;
      else
      {
         if (*iptr <= 1) *kptr = -1;
         else
         {
            ilo = *iptr-1;
            ihi = *iptr;
         }
      }
   }
   if (*kptr == 10)
   {
      if (*xptr < ara[ihi])
      {
         if (*xptr >= ara[ilo]) *kptr = 0;
         else
         {
            istep = 1;
            while (ilo > 1 && *xptr < ara[ilo])
            {
               ihi = ilo;
               ilo = ihi-istep;
               istep *= 2;
            }
            if (ilo <= 1)
            {
               ilo = 1;
               if (*xptr < ara[1]) *kptr = -1;
            }
         }
      }
      else
      {
         istep = 1;
         while (ihi < *iptr && *xptr > ara[ihi])
         {
            ilo = ihi;
            ihi = ilo+istep;
            istep *= 2;
         }
         if (ihi >= *iptr)
         {
            ihi = *iptr;
            if (*xptr > ara[*iptr]) *kptr = 1;
         }
      }
      if (*kptr == 10)
      {
         do
         {
            middle = (ilo+ihi)/2;
            if (*xptr >= ara[middle])
            {
               ilos = ilo;
               ilo = middle;
            }
            else ihi = middle;
         }
         while (middle != ilos);
      }
   }
   if (*kptr == -1) *jptr = 1;
   else
   {
      if (*kptr == 1) *jptr = *iptr;
      else
      {
         *kptr = 0;
         *jptr = ilo;
      }
   }
}

/*------------------------------------------------------------------*/
double zbrent(x1,x2,ind,value)
   int   ind;
   double   x1, x2, value;

   /* using Brent's method, find the root of a function known to lie
    * between x1 and x2; the root returned as x3, will be refined until
    * its accuracy is smaller than 1e-8 times f(x1) or f(x2), depending
    * which one is largest;
    */
{
   double   y1, y2, y, a, b, fa, fb, x3, fc, c, d, e, tol1, xm, s, p,
      q, r, eps;
   eps = 1e-8;
   a = x1;
   b = x2;
   x3 = 0;
   fa = ppvalu(brek,coef,&l,&k,&a,ind)-value;
   fb = ppvalu(brek,coef,&l,&k,&b,ind)-value;
   fc = fb;
   y1 = fabs(a);
   y2 = fabs(b);
   if (y1 > y2) y = y1;
   else y = y2;
   tol1 = 2*eps*y;
   while (x3 == 0)
   {
      if (fb*fc > 0)
      {
         c = a;
         fc = fa;
         d = b-a;
         e = d;
      }
      if (fabs(fc) < fabs(fb))
      {
         a = b;
         b = c;
         c = a;
         fa = fb;
         fb = fc;
         fc = fa;
      }
      xm = .5*(c-b);
      if (fabs(xm) <= tol1 || fb == 0) x3 = b;
      else
      {
         if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
         {
            s = fb/fa;
            if (a == c)
            {
               p = 2*xm*s;
               q = 1-s;
            }
            else
            {
               q = fa/fc;
               r = fb/fc;
               p = s*(2*xm*q*(q-r)-(b-a)*(r-1));
               q = (q-1)*(r-1)*(s-1);
            }
            if (p > 0) q *= -1;
            p = fabs(p);
            y1 = 3*xm*q-fabs(tol1*q);
            y2 = fabs(e*q);
            if (y1 < y2) y = y1;
            else y = y2;
            if (2*p < y)
            {
               e = d;
               d = p/q;
            }
            else
            {
               d = xm;
               e = d;
            }
         }
         else
         {
            d = xm;
            e = d;
         }
         a = b;
         fa = fb;
         if (fabs(d) > tol1) b += d;
         else
         {
            if (xm >= 0) y = tol1;
            else y = -1*tol1;
            b += y;
         }
         fb = ppvalu(brek,coef,&l,&k,&b,ind)-value;
      }
   }
   return(x3);
}
