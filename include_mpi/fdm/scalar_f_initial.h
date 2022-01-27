void ini_rand_field(int * ,double *,double *,double * ,double ,double ,double,ini_power_generator );
void ini_rand_field_2(int * ,double [][3],int *,int ,double,double *,double * ,double ,double ,double,ini_power_generator );
double res_limits(double ,double ,double ,double ,double &,double & );
void initial_hdf5_write_mpi(int *,int *,fdm_psi_mpi ,metric_potential_mpi ,hid_t ,double *,double [][3],double [][3],double ,double,int ,bool );

void initialise_mpi(int * ,int *,fdm_psi_mpi &,metric_potential_mpi &,double [][3],int[] ,double ,double ,double,double,double *,double &,int &,ini_power_generator,gauss_rand_field_gen_mpi,bool);
		   
double ini_power_spec(double );

#include "../../source_mpi/fdm/scalar_f_initial.cpp"
