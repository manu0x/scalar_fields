int evolve_kdk(int *,fdm_psi &,metric_potential &,double [][3],int [],double ,double,double,double,double *,double,int,double dt=1e-4);
int evolve_kdk_openmp(int *,fdm_psi &,metric_potential &,double [][3],int [],double ,double,double,double,double *,double,int,double dt=1e-4);
void evolve_hdf5_write(int *,fdm_psi ,metric_potential ,hid_t ,double ,double,double ,bool );
#include "../source/scalar_f_evolve.cpp"
