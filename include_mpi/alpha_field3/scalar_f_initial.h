
double res_limits(double ,double ,double ,double ,double &,double & );
void initial_hdf5_write_mpi(int *,int *,metric_potential_poisson_mpi ,metric_potential_poisson_mpi,hid_t ,double *,double [][3],double [][3],double ,double,int ,bool );

void initialise_mpi(int * ,int *,metric_potential_poisson_mpi &,field_vel_mpi &,metric_potential_poisson_mpi &,metric_potential_poisson_mpi_ini &,
				double ,int ,double ,double,double,double ,double &,double *,double &,int & ,
								ini_power_generator ,gauss_rand_field_gen_mpi ,bool ,double,double,int,string);
		   
double ini_power_spec(double );

#include "../../source_mpi/alpha_field3/scalar_f_initial.cpp"
