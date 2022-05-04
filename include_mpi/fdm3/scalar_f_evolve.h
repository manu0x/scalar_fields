//int evolve_kdk(int *,int *,fdm_psi_mpi &,metric_potential_mpi &,double [][3],int [],double ,double,double,double,double *,double,int,double ,bool);
int evolve_kdk_openmp(int *,int *,fdm_poisson_mpi &,metric_potential_poisson_mpi &,
				double [][3],int [],double ,double,double,double,double *,double,int,double, bool);
void evolve_hdf5_write(int *,fdm_poisson_mpi,metric_potential_poisson_mpi ,hid_t ,double *,double *,double,double,double,double,int, bool );


#include "../../source_mpi/fdm3/scalar_f_evolve.cpp"
