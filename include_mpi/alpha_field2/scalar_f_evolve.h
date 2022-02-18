int evolve_kdk(int *,int *,field_alpha_mpi &,metric_potential_approx_1_t_mpi &,double [][3],int [],double ,double,double,double,double,double *,double,int,double,int, bool);
int evolve_kdk_openmp(int *,int *,field_alpha_mpi &,metric_potential_poisson_mpi &,double [][3],int [],double ,double,double,double,double,double *,double,int,double,int, bool);


void evolve_hdf5_write(int *,field_alpha_mpi,metric_potential_poisson_mpi ,hid_t ,double *,double,double,double,int,bool );

#include "../../source_mpi/alpha_field2/scalar_f_evolve.cpp"
