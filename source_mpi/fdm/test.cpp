using namespace std;
using namespace std;
#include "../../include_mpi/fdm/include_custom.h"



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
      	int ind[3]{6,3,3};
	
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
		


 	MPI_Type_vector(1,0,n_axis_loc[2]*n_axis_loc[1],MPI_DOUBLE,&c_x_plain);
  	MPI_Type_commit(&c_x_plain);
  
 	//MPI_Type_vector(n_axis_loc[0],n_axis_loc[2],(n_axis_loc[2]+4)*(n_axis_loc[1]+4),MPI_DOUBLE,&c_y_plain);
 	//MPI_Type_commit(&c_y_plain);
  
 	//MPI_Type_vector(n_axis_loc[0]*n_axis_loc[1],1,n_axis_loc[2]+4,MPI_DOUBLE,&c_z_plain);
 	//MPI_Type_commit(&c_z_plain);



///////////////////////////////////////////////////////////////////////////////////////////////////

	int c1=1,c2=1,fail=0;
		
	int tN = n_axis[0]*n_axis[1]*n_axis[2];
	int tN_loc = n_axis_loc[0]*n_axis_loc[1]*n_axis_loc[2];

	int use_omp{1};
	bool use_hdf5_format{true};
	

	
	scalar_field_3d_mpi f(n_axis_loc, cum_lin_ind,false,false);

	int ind_loc[3];
	double fv;
	FILE *fp1,*fp2;
	if(my_corank==0)
	fp1 = fopen("testONE.txt","w");
	else
	fp2 = fopen("testTWO.txt","w");

	int j,k;
	



	 for(i=0;i<n_axis_loc[0];++i)
	 {
		  for(j=0;j<n_axis_loc[1];++j)
		  {
		    for(k=0;k<n_axis_loc[2];++k)
		    {

			ind_loc[0] = i; ind_loc[1] = j; ind_loc[2] = k;
			f.update_field(ind_loc,(double)(i+1));	


		    }

		   }
	}


	for(i=-2;i<n_axis_loc[0]+2;++i)
	 {
		  for(j=0;j<n_axis_loc[1];++j)
		  {
		    for(k=0;k<n_axis_loc[2];++k)
		    {

			ind_loc[0] = i; ind_loc[1] = j; ind_loc[2] = k;
			fv =  f.get_field(ind_loc,give_f);	
			if(my_corank==0)
			fprintf(fp1,"%lf\t",fv);
			else
			fprintf(fp2,"%lf\t",fv);


		    }
			if(my_corank==0)
			fprintf(fp1,"\n");
			else
			fprintf(fp2,"\n");
		   }

		if(my_corank==0)
		fprintf(fp1,"\n\n\n");
		else
		fprintf(fp2,"\n\n\n");
	}
	
	if(my_corank==0)
	fprintf(fp1,"\n###############################################\n");
	else
	fprintf(fp2,"\n###############################################\n");
	
	f.mpi_send_recv();

	for(i=-2;i<n_axis_loc[0]+2;++i)
	 {
		  for(j=0;j<n_axis_loc[1];++j)
		  {
		    for(k=0;k<n_axis_loc[2];++k)
		    {

			ind_loc[0] = i; ind_loc[1] = j; ind_loc[2] = k;
			fv =  f.get_field(ind_loc,give_f);	
			if(my_corank==0)
			fprintf(fp1,"%lf\t",fv);
			else
			fprintf(fp2,"%lf\t",fv);


		    }
			if(my_corank==0)
			fprintf(fp1,"\n");
			else
			fprintf(fp2,"\n");
		   }

		if(my_corank==0)
		fprintf(fp1,"\n\n\n");
		else
		fprintf(fp2,"\n\n\n");
	}

	
	
	
	
	printf("fail is %d\n",fail);
	fftw_mpi_cleanup();
	

	MPI_Finalize();
	
}













