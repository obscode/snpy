extern int spline2(double *, double *, double *, int, int , int , 
   int , int , double , double , int , double ,
   double , int , int , int , int , int , 
   double , double , double , int , int , 
   double , int , int , int , int , 
   int , double **, int *, double **, int *, 
   double *, double *, int *, double *, 
   double *);

extern int evalsp(double *, int , int , int , double *, double *, int ,
            double **);
extern int eval_extrema(int,int,double *,double *,double **, double **, int **,
      int *);
extern int eval_inflect(int , int , double *, double *, double **, double **, double **, 
      int *);
extern double eval_integ(double , double , int , int , double *, double *);
extern int eval_x(double , int , int , double *in_t, double *,double **, int *) ; 
