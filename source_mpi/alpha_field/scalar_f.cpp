using namespace std;
#include "../../include_mpi/alpha_field/include_custom.h"



//////////////////////////////Constants/////////////////////////

///////////////////////////////////////////Todo List///////////////////////////////////////////
//////////////   1) Add error check for successful memory allocation in class scalar_field_3d    /////////////////////////////////////////
//////////////   2) Check cosmo ini conditions as per fdm					 /////////////////////////////////////////
//////////////   3) Check k grid related things							 /////////////////////////////////////////











int main(int argc, char **argv)
{
	clock_t t_start,t_end;

	t_start = clock();


	int mpicheck = 0;
	int num_p = 0;
	int rank;
	mpicheck = MPI_Init(&argc,&argv);
	mpicheck = MPI_Comm_size(MPI_COMM_WORLD,&num_p);
	
	mpicheck = MPI_Comm_rank(MPI_COMM_WORLD,&rank);




	char *param_file_name;
	char fp_sim_info_name[30]("sim_info_"), str_corank[3];

	param_file_name = argv[1];
	printf("Param filename is %s\n",param_file_name);
	

	int i;
	param_alpha p;
	p.load_defaults();

	i =  parser(param_file_name, p);
      	int ind[3]{p.box_n,p.box_n,p.box_n};


	
	fftw_mpi_init();
////////////////////////////////	MPI related...	/////////////////////////////


int nd_cart;
int *my_coords;
int my_corank;
int n_axis[3];
int n_axis_loc[3];
int n_axis_loc_act[3];

int cum_lin_ind;

int * dims,*periods;

MPI_Status stdn,stup;




///////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////	MPI related....	/////////////////////////////////
	n_axis[0]=ind[0];
	n_axis[1]=ind[1];
	n_axis[2]=ind[2];


	

	if(num_p<=8)
		nd_cart = 1;
	else
	{	if(num_p<20)
			nd_cart = 2;
		else
			nd_cart = 3;
	}

	dims = new int [nd_cart];
	periods = new int [nd_cart];
	my_coords = new int [nd_cart];

	for(i=0;i<nd_cart;++i)
	{	periods[i] = 1;
		dims[i]=0;
		


	}

	

	mpicheck = MPI_Dims_create(num_p,nd_cart,dims);

	

	

	mpicheck = MPI_Cart_create(MPI_COMM_WORLD,nd_cart, dims, periods,1,&cart_comm);

	mpicheck = MPI_Cart_get(cart_comm,nd_cart,dims,periods,my_coords);
	
	mpicheck = MPI_Cart_rank(cart_comm,my_coords,&my_corank);

	
	sprintf(str_corank,"%d",my_corank);
	
	strcat(fp_sim_info_name,str_corank);
	strcat(fp_sim_info_name,".txt"); 
	fp_sim_info = fopen(fp_sim_info_name,"w");

	printf("fp info name %s\n",fp_sim_info_name);
	
	fprintf(fp_sim_info,"/////////////////	CART-RANK %d TALKING	//////////////////\n\n",my_corank);
	
	fprintf(fp_sim_info,"num process  %d nd_cart %d dim[0] %d\n",num_p,nd_cart,dims[0]);

	for(i=0;i<3;++i)
	{	
		if(i<nd_cart)
		{
		   double temp_naxis = (((double) n_axis[i])/ ((double) dims[i]));
		   n_axis_loc[i] =  (int)(ceill(temp_naxis));

		
		   

		    if(i==0)
		   {
			
			cum_lin_ind = my_coords[i]*n_axis_loc[i];
			 if(my_coords[i]==(dims[i]-1))
				cum_lin_ind=((dims[i]-1)*n_axis_loc[i]);

			printf("cum_lin_ind %d\n",cum_lin_ind);

	           }
	

		    if((n_axis[i]%dims[i]) != 0)
		     {
		     
			 if(my_coords[i]==(dims[i]-1))
		      		n_axis_loc[i]+=(n_axis[i]-dims[i]*n_axis_loc[i]);

		      }

			
		 

		}

		else
		   n_axis_loc[i] = n_axis[i]; 

		if(num_p==1)
		{
			n_axis_loc[i] = n_axis[i];
			if(i==0)
			cum_lin_ind = 0;

		}
		

		n_axis_loc_act[i] = n_axis_loc[i]+4;

		

	}
	
	fprintf(fp_sim_info,"cO_RANK %d cords %d with nloc 0  %d and cum_ind %d\n",my_coords[0],my_corank,n_axis_loc[0],cum_lin_ind);
		


 	
 	MPI_Type_vector(1,n_axis_loc[2]*n_axis_loc[1],0,MPI_DOUBLE,&c_x_plain);
  	MPI_Type_commit(&c_x_plain);
  
 	//MPI_Type_vector(n_axis_loc[0],n_axis_loc[2],(n_axis_loc[2]+4)*(n_axis_loc[1]+4),MPI_DOUBLE,&c_y_plain);
 	//MPI_Type_commit(&c_y_plain);
  
 	//MPI_Type_vector(n_axis_loc[0]*n_axis_loc[1],1,n_axis_loc[2]+4,MPI_DOUBLE,&c_z_plain);
 	//MPI_Type_commit(&c_z_plain);


///////////////////////////////////////////////////////////////////////////////////////////////////

	int c1=1,c2=1,fail=0;
		
	int tN = n_axis[0]*n_axis[1]*n_axis[2];
	int tN_loc = n_axis_loc[0]*n_axis_loc[1]*n_axis_loc[2];

	field_alpha_mpi f_alpha(n_axis_loc,cum_lin_ind,true,true);
	metric_potential_approx_1_t_mpi phi(n_axis_loc,cum_lin_ind,true,true);
	metric_potential_poisson_mpi poisson_phi(n_axis,n_axis_loc,cum_lin_ind);
	
	//psi.test_ind();
	//psi.test_ind2();
	

	int use_omp{1};
	bool use_hdf5_format{true};
	

	const char name[] = "lcdm_00_pk.dat";
	//printf("\n Building pk...\n");
	ini_power_generator pk(name);
	pk.check(0.0);
	gauss_rand_field_gen_mpi grf(n_axis,n_axis_loc);
	
	//gen.stats_check(3.1);

	

	double a0,ai,omega_dm_ini,Xb_0;
	double k_grid[tN_loc][3],dx[3],dk;
	int kbins,kbin_grid[tN_loc];
	
	set_back_cosmo(a0,ai,Hi,omega_dm_ini,p);
	printf("Hi %lf\nOmega_dm_ini %lf\nai %lf\n h is %lf\n",Hi,omega_dm_ini,ai,h);
	


	initialise_mpi(n_axis,n_axis_loc,f_alpha,phi,poisson_phi,
				k_grid,kbin_grid,a0,ai,Hi,omega_dm_ini,Xb_0,dx,dk,kbins,pk,grf,use_hdf5_format,p.box_length,cum_lin_ind);
	
	
	printf("\nHHHdk is %lf\n",dk);

	

	if(use_omp)
	fail = evolve_kdk_openmp(ind,n_axis_loc,f_alpha,phi,k_grid,kbin_grid,a0,ai,a0,omega_dm_ini,Xb_0,dx,dk,kbins,0.2e-5,cum_lin_ind,use_hdf5_format);
	else
	fail = evolve_kdk(ind,n_axis_loc,f_alpha,phi,k_grid,kbin_grid,a0,ai,a0,omega_dm_ini,Xb_0,dx,dk,kbins,0.2e-4,cum_lin_ind,use_hdf5_format);
	
	
	printf("fail is %d\n",fail);
	fftw_mpi_cleanup();
/**/	
	
	MPI_Finalize();
	
}






