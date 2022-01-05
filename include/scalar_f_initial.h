void ini_rand_field(int * ,double *,double *,double * ,double ,double ,double,ini_power_generator );
void ini_rand_field_2(int * ,double [][3],int *,int ,double,double *,double * ,double ,double ,double,ini_power_generator );
void res_limits(double ,double ,double ,double ,double &,double & );

void initialise(int * ,fdm_psi &,metric_potential &,double [][3],int[] ,double ,double ,double,double,double *,double &,int &,ini_power_generator,gauss_rand_field_gen);
		   
double ini_power_spec(double );

#include "../source/scalar_f_initial.cpp"
