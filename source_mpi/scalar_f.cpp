using namespace std;
#include "../include_mpi/include_custom.h"



//////////////////////////////Constants/////////////////////////

///////////////////////////////////////////Todo List///////////////////////////////////////////
//////////////   1) Add error check for successful memory allocation in class scalar_field_3d    /////////////////////////////////////////
//////////////   2) Check cosmo ini conditions as per fdm					 /////////////////////////////////////////
//////////////   3) Check k grid related things							 /////////////////////////////////////////











int main(int argc, char **argv)
{
	clock_t t_start,t_end;

	t_start = clock();


	int i;
      	int ind[3]{256,256,256};
	
	fftw_mpi_init();
////////////////////////////////	MPI related...	/////////////////////////////

int mpicheck = 0;
int num_p = 0;
int rank;


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


	mpicheck = MPI_Init(&argc,&argv);
	mpicheck = MPI_Comm_size(MPI_COMM_WORLD,&num_p);
	
	mpicheck = MPI_Comm_rank(MPI_COMM_WORLD,&rank);

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

	printf("num proce %d %d %d\n",num_p,nd_cart,dims[0]);

	mpicheck = MPI_Cart_create(MPI_COMM_WORLD,nd_cart, dims, periods,1,&cart_comm);

	mpicheck = MPI_Cart_get(cart_comm,nd_cart,dims,periods,my_coords);
	
	mpicheck = MPI_Cart_rank(cart_comm,my_coords,&my_corank);


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
	
		printf("cO_RANK %d cords %d with nloc 0  %d and cum_ind %d\n",my_coords[i],my_corank,n_axis_loc[0],cum_lin_ind);
		


 	MPI_Type_vector(n_axis_loc[1],n_axis_loc[2],n_axis_loc[2]+4,MPI_DOUBLE,&c_x_plain);
  	MPI_Type_commit(&c_x_plain);
  
 	MPI_Type_vector(n_axis_loc[0],n_axis_loc[2],(n_axis_loc[2]+4)*(n_axis_loc[1]+4),MPI_DOUBLE,&c_y_plain);
 	MPI_Type_commit(&c_y_plain);
  
 	MPI_Type_vector(n_axis_loc[0]*n_axis_loc[1],1,n_axis_loc[2]+4,MPI_DOUBLE,&c_z_plain);
 	MPI_Type_commit(&c_z_plain);



///////////////////////////////////////////////////////////////////////////////////////////////////

	int c1=1,c2=1,fail=0;
		
	int tN = n_axis[0]*n_axis[1]*n_axis[2];
	int tN_loc = n_axis_loc[0]*n_axis_loc[1]*n_axis_loc[2];

	fdm_psi_mpi psi(n_axis_loc,cum_lin_ind,true);
	metric_potential_mpi phi(n_axis,n_axis_loc,true);

	int use_omp{1};
	bool use_hdf5_format{true};
	

	const char name[] = "lcdm_00_pk.dat";
	printf("\n Building pk...\n");
	ini_power_generator pk(name);
	pk.check(0.0);
	gauss_rand_field_gen_mpi grf(n_axis,n_axis_loc);
	
	//gen.stats_check(3.1);

	

	double a0,ai,omega_dm_ini;
	double k_grid[tN_loc][3],dx[3],dk;
	int kbins,kbin_grid[tN_loc];
	
	set_back_cosmo(a0,ai,Hi,omega_dm_ini);
	printf("Hi %lf\nOmega_dm_ini %lf\nai %lf\n",Hi,omega_dm_ini,ai);

	initialise_mpi(n_axis,n_axis_loc,psi,phi,k_grid,kbin_grid,a0,ai,Hi,omega_dm_ini,dx,dk,kbins,pk,grf,use_hdf5_format,cum_lin_ind);
	psi.mpi_send_recv();
	
	//printf("\ndk is %lf\n",dk);
	
/*	if(use_omp)
	fail = evolve_kdk_openmp(ind,psi,phi,k_grid,kbin_grid,a0,ai,a0,omega_dm_ini,dx,dk,kbins,0.4e-4,use_hdf5_format);
	else
	fail = evolve_kdk(ind,psi,phi,k_grid,kbin_grid,a0,ai,a0,omega_dm_ini,dx,dk,kbins,0.4e-4,use_hdf5_format);
	
	printf("fail is %d\n",fail);

	fftw_mpi_cleanup();
*/	

	MPI_Finalize();
	
}






