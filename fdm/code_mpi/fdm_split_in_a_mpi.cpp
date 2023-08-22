

#include <stdio.h>
#include <math.h>

#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <hdf5.h>


#include <limits>


#include <stdlib.h>
#include <string.h>
#include <string>

#include <algorithm>
#include "../../imex/imex_classes.cpp"
#include "../../GPE/GPE_classes_split_exp_uni_mpi.cpp"
#include "./utilities_mpi.cpp"

using namespace std;

#define pie M_PI



//////////GLobal constants/////////////



///////////////////////////////////////








void initialise_kgrid(double *k,double dx,int N)
{

	int i,j,l,ci; double di,dj,dl,dci,dk;
	dk = 2.0*pie/(((double)(N))*dx);
	
	for(i=0;i<N;++i)
	{
	   di = (double)i;

	   if(i<(N/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) N))*dk;


	}






}


double run(double da,int dim,double *mass_err,int calvel,int cal_da_cons,int prntfp,int prnt,int argc,char **argv)
{

	int my_rank;
	int mpicheck;
	int N;
	int myNx,myN_tot,N_tot,cum_lin_ind;

	mpicheck = MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	int a_steps;
	double box_len,a,a_end,a_start,ak,da_up,xval,dbi,dbj,dbk,dx;
	double HbyH0;
	double x0[3],xv[3];

	double velmax;

	da_up=da;

	//fftw_init_threads();


/////////////////////////////// Box & res. setting ///////////////////////////////////////////////
	box_len = 1.0;
	x0[0]=-0.5*box_len; x0[1]=-0.5*box_len; x0[2]=-0.5*box_len;
	
	//n  = 2.0/box_len;
	//N = 256;
	N = atoi(argv[2]);
	dx = box_len/(double(N));

	int N2,N3;
	N2 = N*N;
	N3=N*N*N;
	
	//N = ((int)(box_len/dx));

	printf("N %d\n",N);

////////////	Imex declare and data write file declare ////////////////////////


	

////////////////////////	Declare field /////////////////////////////////////////////
char *f1paramfile;
char *fpsi_ini_file;

char pas[20] = {"psi"};

GPE_field_mpi  psi_1(dim,N,0,my_rank,pas,1,16); 
myNx = psi_1.myNx;
myN_tot=psi_1.myN_tot;
cum_lin_ind = psi_1.cum_lin_ind;

f1paramfile = argv[3];
fpsi_ini_file = argv[4];

psi_1.read_param_from_file(f1paramfile);
	

psi_1.print_params_set_kappa();







////////////////////////////// Time & da settings ////////////////////////////////////////////////

	a_start = psi_1.ai;
	a_end = 1.0;
	
	a_steps = (int)((a_end-a_start)/da);
	
	if(prnt&&(my_rank==0))
	printf("da %lf N %d\n",da,N);
	HbyH0= psi_1.HbyH0(a_start);
	if(my_rank==0)
	printf("chking c %lf %lf\n",(a_start*a_start*HbyH0),psi_1.kppa/(a_start*a_start*a_start*HbyH0));


	
/////////////////////////////////////////RK things/////////////////////////////////////////
	int i,j;

////////////////////////////// k things /////////////

	double k_grid[N],ksqr;


//////////////////////////////  variables  /////////////////////////////////////////////////////



	double lambda;



    int ind_loc[3] = {myNx,N,N};
	

	initialise_kgrid(k_grid,dx,N);

	//psi_1.initialise_random(k_grid) ; 
	//psi_1.read_from_initial();

	read_psi_from_hdf5_mpi(fpsi_ini_file,psi_1.psi,ind_loc, cum_lin_ind);
	printf("INI file is %s\n",fpsi_ini_file);
	
	//printf("no of threads %d\n",omp_get_num_threads());
	//omp_set_num_threads(4);
	//printf("no of threads %d\n",omp_get_num_threads());

	
	

	if(prnt)
	printf("Initialization done from rank %d\n",(my_rank));
	
	psi_1.set_field();

	//psi_1.write_hdf5_mpi(1.0/a_start -1.0);
	


///////////////////////  Evolution ///////////////////////////////////////////////////////////////
    int ii,jj,kk;
	int s_cntr,acntr,printcntr,fpcntr,err_cntr=1,fail=0;

	ii=-1;
	jj=-1;
	kk=0;


	printcntr = (int)(((double) a_steps)/100.0);
	fpcntr = (int)(((double) a_steps)/10.0);

	if(a_steps<=100)
	{	printcntr = 1.0;
	
	}	//printf("%d\n",printcntr);

	//Variables for holding temporary values
	double c_psi[2],c_psi_amp2,delta;
	
	double delta_sum,fac,Vmax=-0.00000000000001,Vmax_mpi_glbl;
	double da_constraints[2];
	

	
	if(my_rank==0)
	printf("Starting Run..,\n");

	double cor_fac;
	
	for(a=a_start,acntr=0;(a<=a_end)&&(!fail)&&(1);a+=da,++acntr)
	{




		///////////////	Data Writing ////////////////////////////////

		if((acntr%fpcntr==0)&&(prntfp))
			psi_1.write_hdf5_mpi(1.0/a -1.0);


			
		for ( i = 0; i <myN_tot ; i++)
		{
			if(Vmax<fabs(psi_1.V_phi[i][0]))
			{			Vmax = fabs(psi_1.V_phi[i][0]);
				
			}
		}

		if(cal_da_cons)
		{ 
		  MPI_Allreduce(&Vmax, &Vmax_mpi_glbl, 1, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);	
		  da_constraints[0] = a*a*psi_1.HbyH0(a)*2.0*M_PI/Vmax_mpi_glbl;
		  da_constraints[1] = 4.0*a*a*a*psi_1.HbyH0(a)*(dx*dx)/(3.0*M_PI*psi_1.kppa);

		}
	
		if((acntr%printcntr==0)&&prnt)
		{	if(calvel)
			velmax = psi_1.calc_v(k_grid);
			printf("time %lf %lf %lf  da %lf  from rank %d\n",a/a_end,psi_1.mass,velmax,da,my_rank);

			if(cal_da_cons)
			printf("Constraints da<  {%lf,%lf}\n",da_constraints[0],da_constraints[1]);
			if(calvel)
			printf("Constraints on (%lf)dx <%lf\n",dx,M_PI/velmax);

			//printf("SV da<=  %lf\n",a*a*a*psi_1.HbyH0(a)*(dx*dx)*(4.0/(3.0*pie))/psi_1.kppa);

		}


	  psi_1.split_step(a,da,k_grid,acntr);
		
	}///// ENd f Time Evolution /////////////////////////

	//fclose(fp);
	//fftw_cleanup_threads();

	*mass_err = psi_1.max_mass_err;
	*(mass_err+1)= psi_1.max_mass_err;


	if(prnt)
	printf("N %d\n Run en los %lf abs err %lf from rank %d\n",N,*(mass_err),*(mass_err+1),my_rank);

	psi_1.write_hdf5_mpi(1.0/a -1.0);

	
	return(psi_1.max_mass_err);








}






int main(int argc, char ** argv)
{


	clock_t t_start,t_end;
	double start_t_mpi,end_t_mpi;

	t_start = clock();


	int mpicheck = 0;
	int num_p = 0;
	int rank;
	mpicheck = MPI_Init(&argc,&argv);
	mpicheck = MPI_Comm_size(MPI_COMM_WORLD,&num_p);
	
	mpicheck = MPI_Comm_rank(MPI_COMM_WORLD,&rank);


	MPI_Barrier(MPI_COMM_WORLD);
	start_t_mpi = MPI_Wtime();

	//fftw_mpi_init();


	double da = 3e-4;
	double dx = 4e-3;
	double mass_err[2],mass_loss;

	double dx_l=2e-3,dx_u = 4e-2;
	double da_l= 1.25e-4,da_u = 1e-3;

	//double ddx = (dx_u-dx_l)/(20.0);
	//double dda = (da_u-da_l)/(20.0);

	FILE *fp;
	if(rank==0)
	fp = fopen("imex_ft_mpi.txt","w");
	

	//for(da=da_l;da<=da_u;da*=2.0)
	{


		//for(dx = dx_l;dx<=dx_u;dx+=ddx)
		{
			//dx = 2e-2;
			da = atof(argv[1]);
		
			//printf("\n\n RUNNING da = %lf\n\n",da);
		/*	if(da==da_l)
			mass_loss = run(da,3,512,mass_err,1,0,argc,argv); //(double da,int dim,int N,double *mass_err,int prntfp,int prnt,int argc,char **argv)
			else
		*/	mass_loss = run(da,3,mass_err,0,1,1,1,argc,argv);
			printf("%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,da,da/(dx*dx),mass_loss,*mass_err,*(mass_err+1));
			if(rank==0)
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,da,da/(dx*dx),mass_loss,*mass_err,*(mass_err+1));
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
	end_t_mpi = MPI_Wtime();
	t_end = clock();

	printf("\nMPI time from rank %d is: %lf\n",rank,end_t_mpi-start_t_mpi);
	printf("\nClock time from rank %d is: %lf\n",rank,t_end-t_start);

	if(rank==0)
	fclose(fp);

	mpicheck = MPI_Finalize();




}
