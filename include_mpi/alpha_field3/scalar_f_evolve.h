/*int evolve_kdk(int *,int *,field_alpha_mpi &,metric_potential_approx_1_t_mpi &,double [][3],int [],double ,double,double,double,double,double *,double,int,double,int, bool);*/
int evolve_kdk_openmp(int *,int *,metric_potential_poisson_mpi &,field_vel_mpi &,metric_potential_poisson_mpi &,double [][3],int [],double ,double,double,double,double,double *,double,int,double,int, bool);


void evolve_hdf5_write(int *,metric_potential_poisson_mpi ,metric_potential_poisson_mpi ,hid_t ,double *,double ,double ,double ,double ,
														int ,bool );

#include "../../source_mpi/alpha_field3/scalar_f_evolve.cpp"
