# include <stdio.h> 
# include <math.h> 
# include <stdlib.h> 
#include <string.h> 

#define MXPTS 200

double bisec(double , double);
double *center_dm15(double , double * , int) ;
double weight[200];
double Wdm15(char *, double, double, double); 
double triweight(double, double, double,  double);
double trifunction(double, double, double, double);
double logfunction(double, double, double, double, double);
double width(double);  
double *zeroes(double *,int); 
int dm15temp(double, int, double *, double *, double *, double *, 
      double *, double *, double *, double *, double *, double *, 
      double *, double *, double *, int *, int, char *);


struct arrays {
  char temp[40];
  char sigma[40];
} names[40];

/*
  Program to construct an interpolated type Ia light curve template
  (in dm15) from the original templates in the 'template set'.
  Information of the templates is in the file 'dm15ṫemps.data'. The
  weighting function used is a triangle function
  
  PRS, March 2006
*/

int main(int argc, char **argv) {
   double t[MXPTS], tsigma[MXPTS], BB[MXPTS], VV[MXPTS], RR[MXPTS], II[MXPTS],
      BBsig[MXPTS], VVsig[MXPTS], RRsig[MXPTS], IIsig[MXPTS];
   double shiftV, shiftR, shiftI;
   char path[] = ".";
   int npts, i, method;
   double dm15;
  if(argc < 2) {
    printf("Usage:  dm15temp   dm15 [method]\n\n");
    exit(0);
  }
  
  /* dm15 wanted for the output template */
  dm15 = atof(argv[1]);
  if ( argc > 2) {
     method = atoi(argv[2]);
  } else {
     method = 1;
  }
  
  /* Allowed range of dm15 */
  if(dm15 > 1.93 || dm15 < 0.83) {
    printf("\nAllowed range: 1.93 > dm15 > 0.83\n\n");
    exit(0);
  }

  dm15temp(dm15, MXPTS, t, tsigma, BB, BBsig, VV, VVsig, RR, RRsig, 
        II, IIsig, &shiftV, &shiftR, &shiftI, &npts, method, path);

  /* Printing out results */ 
  printf("# Constructed light curve template for dm15=%.3f mag\n",dm15);  
  printf("# Columns: \n#  1    t(Bmax) \n#  2-3  B-B(max)   sigma[B-B(max)]^2  \n#  4-5  V-V(max)   sigma[V-V(max)]^2 \n#  6-7  R-R(max)   sigma[R-R(max)]^2 \n#  8-9  I-I(max)   sigma[I-I(max)]^2 \n");
  
  for(i=0; i<npts; i++) {
    printf("%5.1lf %10.3lf %10.6lf %10.3lf %10.6lf %10.3lf %10.6lf %10.3lf %10.6lf\n",t[i],BB[i],BBsig[i],VV[i]-shiftV,VVsig[i],RR[i]-shiftR,RRsig[i],II[i]-shiftI,IIsig[i]); 
  }
}


int dm15temp(double dm15, int size, double *tt, double *ttsig, double *BB,
      double *BBsig, double *VV, double *VVsig, double *RR, double *RRsig,
      double *II, double *IIsig, double *shiftV, double *shiftR, 
      double *shiftI, int *npts, int method, char *path) { 
  int i=0, j=0, k=0, N, JMAX=45, M, l, tempnum;
  /* character buffer for the fully qualified file name */
  char filename[1024];
  double offsetV[20],offsetR[20],offsetI[20],dm15_templates[20];
  double temp3,temp4,temp5,temp6,temp7;
  double center_min, center_max;
  double f1,f2,fmid1,fmid2, W[MXPTS], rtb, dx, xmid, *w1, ddm15, ddm15min;
  float *w2;
  double t[MXPTS], tsigma[MXPTS], B[MXPTS], V[MXPTS], R[MXPTS], I[MXPTS], 
         Bsigma[MXPTS], Vsigma[MXPTS], Rsigma[MXPTS], Isigma[MXPTS];
  char *temp1[20], temp2, outfile1[20], outfile2[20];
  FILE *arch1, *arch3, *arch4, *arch5;
  FILE *arch;
  int nchars, imin;
  
  *shiftV=0; *shiftR=0; *shiftI=0;

  /* Reading data of the templates */
  nchars = snprintf(filename, 1024, "%s/%s", path, "dm15temps.dat"); 
  if (nchars > 1024) {
     fprintf(stderr, "Error:  path name for templates file too long\n");
     return(-1);
  }
  arch1 = fopen(filename,"r");
  while(fscanf(arch1,"%s %s %lf %lf %lf %lf", names[i].temp, names[i].sigma,
        &temp3, &temp4, &temp5, &temp6)!=EOF) {
    dm15_templates[i] = temp3;
    offsetV[i] = temp4;
    offsetR[i] = temp5;
    offsetI[i] = temp6;
    i++;
  }
  N=i; i=0; j=0;
  fclose(arch1);

  if (method == 0) {
     /* Just take the closest match */
     imin = -1;
     for (i = 0 ; i < N ; ++i) {
        ddm15 = fabs(dm15 - dm15_templates[i]);
        if (imin == -1) {
           ddm15min = ddm15;
           imin = 0;
        } else {
           if (ddm15 < ddm15min) {
              ddm15min = ddm15;
              imin = i;
           }
        }
     }
     for (i = 0 ; i < N ; ++i ) {
        if (i!=imin) { 
           weight[i] = 0.0;
        } else {
           weight[i] = 1.0;
        }
     }
     w1 = weight;
     for (i = 0 ; i < N ; ++i) {
        fprintf(stderr, "w(%f) = %f\n", dm15_templates[i], w1[i]);
     }
  } else {
    if(dm15 < 0.90) {
      center_min= dm15 - 0.10;
      center_max= dm15 + 0.05;
    }
    else {
      center_min= dm15 - 0.15;
      center_max= dm15 + 0.20;
    }
 
    /*  
        Bisection method to find the weights of a given dm15 template
        changing the central value of the triangle function
    */
 
    i=0;
    w1 = center_dm15(center_min, dm15_templates, N);
    f1= w1[N];
    w1 = center_dm15(center_max, dm15_templates, N);
    f2= w1[N];
  
    fmid1=bisec(f1,dm15);
    fmid2=bisec(f2,dm15);
    
    rtb=fmid1 < 0.0 ? (dx=center_max-center_min,center_min): (dx=center_min-
                       center_max,center_max);
    for(j=1; j <= JMAX; j++) {  
      xmid = rtb + (dx *= 0.5);
      /* w1 is the vector of weights */
      w1=center_dm15(xmid,dm15_templates,N);
      fmid1=w1[N];
      fmid2=bisec(fmid1,dm15);
      if(fmid2 <= 0.0)  rtb=xmid;   
      if((fabs(dx) < 0.0000001) || (fmid2 == 0.0)) {j=JMAX+1;}
    }
  }

  /* 
     Finding the difference between the magnitude at maximum and the
     magnitude at Bmax for the new template.
  */

  for(j=0; j < N; j++) {
    W[j]=w1[j];
    *shiftV = offsetV[j] * W[j] + *shiftV;
    *shiftR = offsetR[j] * W[j] + *shiftR;
    *shiftI = offsetI[j] * W[j] + *shiftI;
  } 
  for(j = 0 ; j < size ; ++j) {
     tt[j] = ttsig[j] = 0;
     BB[j] = BBsig[j] = VV[j] = VVsig[j] = 0;
     RR[j] = RRsig[j] = II[j] = IIsig[j] = 0;
  }
  j=0; i=0; k=0;
  
  /* 
     Weighting the templates in the template set 
  */
  
  for(l=0; l < N; l++) {
    if(W[l] > 0) {
      nchars = snprintf(filename, 1024, "%s/%s", path, names[l].temp);
      if (nchars > 1024){
         fprintf(stderr, "Error: filename for template too long\n");
         return(-1);
      }
      arch3=fopen(filename,"r");
      
      /* Read template data */ 

      while(fscanf(arch3,"%lf %lf %lf %lf %lf",&temp3,&temp4,&temp5,&temp6,
            &temp7)!=EOF){
        t[j]=temp3; tt[j]=temp3;
        B[j]=temp4; 
        V[j]=temp5; 
        R[j]=temp6;
        I[j]=temp7;
        j++;
      }
      M=j;

      if (M > size) M = size;
      *npts = M;
      fclose(arch3);

      
      nchars = snprintf(filename, 1024, "%s/%s", path, names[l].sigma);
      if (nchars > 1024){
         fprintf(stderr, "Error: filename for template too long\n");
         return(-1);
      }
      arch4=fopen(filename,"r");

      /* Read the approx. uncertainties in each template (saved as look-up tables)  */

      j = 0;
      while(fscanf(arch4,"%lf %lf %lf %lf %lf",&temp3, &temp4, &temp5, &temp6, &temp7)!=EOF){
        tsigma[j]=temp3;  ttsig[j] = temp3;
        Bsigma[j]=temp4;
        Vsigma[j]=temp5;
        Rsigma[j]=temp6;
        Isigma[j]=temp7;
        j++;
      }

      j=0;
      fclose(arch4);

      /* Weighting */

      for(i=0; i<M; i++) {
        BB[i] = B[i] * W[l] + BB[i]; 
        VV[i] = V[i] * W[l] + VV[i]; 
        RR[i] = R[i] * W[l] + RR[i];
        II[i] = I[i] * W[l] + II[i];
        BBsig[i] = pow((Bsigma[i] * W[l]),2.) + BBsig[i]; 
        VVsig[i] = pow((Vsigma[i] * W[l]),2.) + VVsig[i]; 
        RRsig[i] = pow((Rsigma[i] * W[l]),2.) + RRsig[i];
        IIsig[i] = pow((Isigma[i] * W[l]),2.) + IIsig[i];         
      }          
    }    
  }
  return(0);
}

/* Difference function  */
double bisec(double center,double dm15) {
  double f;
  f = center - dm15;
  return f;
}


/* This function returns the final vector of weights */

double *center_dm15(double center, double dm15_templates[],int N) {
  int j=0,i=0;  
  double WE=0, sum=0, DM15=0, *wit, functionwidth=0, dm15_temp_width=0, dm15_temp=0;
  
  for(i=0; i<200;i++){
    weight[i]=0;
  }
  i=0;
  for( j=0; j < N; j++) {  
    dm15_temp = dm15_templates[j];
    functionwidth = width(center);
    WE = Wdm15("triangle",dm15_temp, center, functionwidth);
    weight[j] = WE;
    sum = sum + WE;
  }

  for(i=0; i<N; i++){
    weight[i] = (weight[i]/sum);
  }
  for(i=0; i<N; i++) {
    DM15 = DM15 + weight[i]*dm15_templates[i];
  }
  weight[N]= DM15;
  return weight;
}

/*
  Definition of different weighting functions: triangle, box function,
  and a log function (the log function doesn't work !)
*/

double Wdm15(char FUNC[], double center_temp, double center_dm15, double fwidth) { 
  int array[3],i=0,j,C;
  double weight;
  struct functions
  {
    char func1[50];
    char func2[50];
    char func3[50];
  };
  struct functions name={"triangle","box","log"};
  C = strcmp(FUNC,name.func1);
  array[0]=C;
  C = strcmp(FUNC,name.func2);
  array[1]=C;
  C = strcmp(FUNC,name.func3);
  array[2]=C;
  for(i=0; i < 3; i++) {
    if(array[i] == 0) {j=i+1;}
  }
  weight=triweight(center_temp, center_dm15, fwidth, j);
  return weight;
}

/* 
   Defining the weights under the different weighting functions.
   Only the triangle function is implemented now
*/

double triweight(double dm15_temp, double center,  double w, double function) {
  double fmin,fmax,wt;
  fmin = center - w/2;  
  fmax = center + w/2;
  if ( (dm15_temp < fmin) || (dm15_temp > fmax) ){
    return 0;
  }
  if( dm15_temp <= center ){
    wt = trifunction(dm15_temp,w,center,0);
    return wt;
  }
  if( dm15_temp > center ){
    wt = trifunction(dm15_temp, w, center, 1);
    return wt;
  }

}

/* Triangle function */ 

double trifunction(double x, double w, double c, double sw) { 
  double triangle;  
  if(sw == 0) {
    triangle = (2/w)*(x-c) + 1;
    return triangle;
  }
  if(sw == 1) {
    triangle = -(2/w)*(x-c) +1;
    return triangle;
  }
}

/* Log function */
double logfunction(double xmin, double xmax, double center, double w, double sw) { 
  double integral;  
  if(sw == 0) {
         integral = - (xmax - xmin)*(1/(center-w/2)) + log(xmax/xmin);
        return integral;
    }
    if(sw == 1) {
         integral = - (xmax - xmin)*(1/(center+w/2)) + log(xmax/xmin);
        return integral;
    }
}

/* 
   The different widths of the triangle function in the dm15 axis are defined. 
   The widths are user defined, depending on the template set, and the 
   coverage (sampling) of dm15 values they have  
*/

double width(double center) {
  double delta1, delta3, delta2, f1 , f2, f3, res;
  delta1=0.20;
  delta3=0.35;
  delta2=0;  
  f1=1.0;
  f2=1.4;
  f3=0;
  if(center <= f1) {
    res=delta1;
    return res;
  }
  if( (center > f1) && (center <= f2) ) { 
    delta2 = ((delta3 - delta1)/(f2-f1))*(center - f1) + delta1;
    res = delta2;  
    return res;
  }
  if(center > f2) {
    res=delta3;
    return res;
  }
}
 
double *zeroes(double vector[],int N) {
  int i=0;
  for(i=0; i<N; i++){
    vector[i]=0;
  }
  return vector;
}
