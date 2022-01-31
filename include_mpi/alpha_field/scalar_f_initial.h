
double res_limits(double ,double ,double ,double ,double &,double & );
void initial_hdf5_write_mpi(int *,int *,field_alpha_mpi ,metric_potential_approx_1_t_mpi,hid_t ,double *,double [][3],double [][3],double ,double,int ,bool );

void initialise_mpi(int * ,int *,field_alpha_mpi &,metric_potential_approx_1_t_mpi &,metric_potential_poisson_mpi &,
				double ,int ,double ,double,double,double ,double *,double &,int & ,
								ini_power_generator ,gauss_rand_field_gen_mpi ,bool ,double,int );
		   
double ini_power_spec(double );

#include "../../source_mpi/alpha_field/scalar_f_initial.cpp"
