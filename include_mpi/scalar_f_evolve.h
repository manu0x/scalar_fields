int evolve_kdk(int *,fdm_psi_mpi &,metric_potential_mpi &,double [][3],int [],double ,double,double,double,double *,double,int,double ,bool);
int evolve_kdk_openmp(int *,fdm_psi_mpi &,metric_potential_mpi &,double [][3],int [],double ,double,double,double,double *,double,int,double, bool);
void evolve_hdf5_write(int *,fdm_psi_mpi ,metric_potential_mpi ,hid_t ,double *,double,double,bool );
#include "../source_mpi/scalar_f_evolve.cpp"
