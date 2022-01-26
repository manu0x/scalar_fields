int evolve_kdk(int *,int *,field_alpha_mpi &,metric_potential_approx_1_t_mp &,double [][3],int [],double ,double,double,double,double *,double,int,double ,bool);
int evolve_kdk_openmp(int *,int *,field_alpha_mpi &,metric_potential_approx_1_t_mp &,double [][3],int [],double ,double,double,double,double *,double,int,double, bool);
void evolve_hdf5_write(int *,field_alpha_mpi,metric_potential_approx_1_t_mpi ,hid_t ,double *,double,double,bool );

#include "../../source_mpi/alpha_field/scalar_f_evolve.cpp"
