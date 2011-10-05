#define n_filters 16
#define max_points 5000
extern double x[n_filters][max_points], y[n_filters][max_points],
             z[n_filters][max_points], wz[n_filters][max_points];
extern char name[n_filters][max_points][80];
extern int np[n_filters];
int load_data(char *);
int dm15temp(int, double, double *, int, double *, double *, int, double, double,
      double, double);
