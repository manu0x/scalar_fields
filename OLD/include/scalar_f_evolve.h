int evolve_kdk(int *,fdm_psi &,metric_potential &,double [][3],int [],double ,double,double,double,double *,double,int,double ,bool);
int evolve_kdk_openmp(int *,fdm_psi &,metric_potential &,double [][3],int [],double ,double,double,double,double *,double,int,double, bool);
void evolve_hdf5_write(int *,fdm_psi ,metric_potential ,hid_t ,double *,double,double,bool );
#include "../source/scalar_f_evolve.cpp"
