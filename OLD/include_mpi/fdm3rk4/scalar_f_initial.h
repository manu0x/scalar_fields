void ini_rand_field(int * ,double *,double *,double * ,double ,double ,double,ini_power_generator );
void ini_rand_field_2(int * ,double [][3],int *,int ,double,double *,double * ,double ,double ,double,ini_power_generator );
double res_limits(double ,double ,double ,double ,double &,double & );
void initial_hdf5_write_mpi(int *,int *,fdm_poisson_mpi ,metric_potential_poisson_mpi  ,hid_t ,double *,double **,double **,double ,double,int ,bool );

void initialise_mpi(int * ,int *,fdm_poisson_mpi &,metric_potential_poisson_mpi &,metric_potential_poisson_mpi_ini &,
		double **,int* ,double ,double ,double,double,double *,double &,int &,ini_power_generator,gauss_rand_field_gen_mpi,bool,double,int);
		   
double ini_power_spec(double );

#include "../../source_mpi/fdm3rk4/scalar_f_initial.cpp"
