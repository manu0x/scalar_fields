

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <mpi.h>
#include <limits>


#include <stdlib.h>
#include <string.h>
#include <string>

#include <algorithm>
#include "../imex/imex_classes.cpp"
#include "../GPE/GPE_classes.cpp"

using namespace std;

#define pie M_PI



//////////GLobal constants/////////////



///////////////////////////////////////








void initialise(double *k,double dx,int N)
{

	int i,j,l,ci; double di,dj,dl,dci,dk;
	dk = 2.0*pie/(((double)(N))*dx);
	
	for(i=0;i<N;++i)
	{
	   di = (double)i;

	   if(i<=(N/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) N))*dk;


	}






}




double run(double dt,int N,double *mass_err,int argc,char **argv,int prntfp,int prnt)
{

	int t_steps;
	double box_len,t_end,t_start,xval,dbi,dbj,dbk,dx;
	double xv,xini;

	fftw_init_threads();
	
///////////////////////////////File for data/////////////////////////////////////////////////////
/*
	FILE *fp = fopen("data_ft.txt","w");
	FILE *fp2 = fopen("data2_ft.txt","w");
	FILE *fplap = fopen("lap_ft.txt","w");
	FILE *fpmass = fopen("mass_ft.txt","w");
	FILE *fptime = fopen("tm_ft.txt","w");
*/
/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	


   




/////////////////////////////// Box & res. setting ///////////////////////////////////////////////
	
	xini= -20.0;
	box_len = 40.0;
	//n  = 2.0/box_len;
	N = 256;
	dx = box_len/(double(N));
	
	//N = ((int)(box_len/dx));
	//int N3 = N*N*N;
	
	printf("N %d\n",N);
////////////////////////////// Time & dt settings ////////////////////////////////////////////////

	t_start = 0.0;
	t_end = 5.0;
	
	t_steps = (int)((t_end-t_start)/dt);
	
	if(prnt)
	printf("dt %lf N %d\n",dt,N);


	
/////////////////////////////////////////RK things/////////////////////////////////////////
	int i,j;




////////////////////////////// Psi variables  /////////////////////////////////////////////////////



	double lambda;


	char fp_name[30]("data_1d_");

	char fp_phi1_r[20]("phi1_r_ini.txt");
	char fp_phi1_i[20]("phi1_i_ini.txt");

	
	
	
    int stages;
	char *imex_file; 
    char *f1paramfile;
	char *f2paramfile;

    stages = atoi(argv[1]);
	imex_file = argv[2];
    f1paramfile = argv[3];
	f2paramfile = argv[4];
	
	strcat(fp_name,imex_file);
	printf("ImEx table filename is %s and stages given by u is %d and out file %s \n",imex_file,stages,fp_name);

	FILE *fp = fopen(fp_name,"w");
	

	imex_table imx(stages);
    
 
  	imx.read_from_file(imex_file);
  	imx.print_table();


    GPE_field_1d  psi_1(0,N,imx.s,8);
  

    psi_1.read_from_file(f1paramfile);


    psi_1.print_params();
	

	psi_1.initialise(dx,xini);
	
	
	
	
	double k_grid[N];

	initialise(k_grid,dx,N);

	if(prnt)
	printf("Initialization done\n");
	
	psi_1.set_field();


///////////////////////  Evolution ///////////////////////////////////////////////////////////////
    
	int s_cntr,tcntr,printcntr,fpcntr,err_cntr=1,fail=0;

	printcntr = (int)(((double) t_steps)/100.0);
	fpcntr = (int)(((double) t_steps)/100.0);
	if(t_steps<=100)
	{	printcntr = 1.0;
	
	}	//printf("%d\n",printcntr);
	double t,vel_val[2],c_psi[2],c_psi_amp2,Vval,amp,avg_amp;
	
	double fdt,amp_ini;
	

	double kv;
	int ind;
	double sol[2],theta,fs;

	
	

	psi_1.reset_consv_quant(1);
	
	printf("Starting Run..,\n");

	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(10);t+=dt,++tcntr)
	{

		
		if((tcntr%printcntr==0)&&prnt)
		printf("time %lf  mass %.10lf\n",t/t_end,psi_1.mass);

		psi_1.do_forward_fft();
		

		psi_1.reset_consv_quant();
		

	 	for(i=0;i<(N);++i)
		{	

			

			dbi = (double)(i);
			xv = xini + dx*dbi;
			if(tcntr%err_cntr==0)
			{
				psi_1.cal_err_at_point(i,xini,dx,t,!tcntr);
				

			}
			
			if((tcntr%fpcntr==0)&&(prntfp))
			{
				dbi = (double)(i);
			   xv = xini+dbi*dx;
			   theta = 0.5*psi_1.c*(xv-psi_1.x0)-(0.25*psi_1.c*psi_1.c -psi_1.a)*t;
	   		   fs = (2.0*psi_1.a/psi_1.bta)/cosh(sqrt(psi_1.a)*(xv-psi_1.x0-psi_1.c*t));
			   sol[0] = cos(theta)*fs;
			   sol[1] = sin(theta)*fs;
				
				fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",xv,psi_1.psi[i][0],psi_1.psi[i][1],sol[0],sol[1]);


			}
	
			lambda =(k_grid[i]*k_grid[i])*imx.im_a[0][0]*dt;
			//printf("lmda  %lf\n",lambda);
			kv = k_grid[i];		

			psi_1.update_fft_fields(i,kv,lambda);
			

		}

		if(tcntr%err_cntr==0)
		{	psi_1.conserve_err();
			

		}
		
		if((tcntr%fpcntr==0)&&(prntfp))
		fprintf(fp,"\n\n\n");
		
		//printf("chk %lf %lf %lf %lf\n",psi_1.psi[100][1],psi_1.fpGpsi[100][1],psi_1.psi[1100][1],psi_1.fpGpsi[1100][1]);
		psi_1.do_back_fft();
		
		//printf("chk %lf %lf %lf %lf\n",psi_1.psi[100][1],psi_1.fpGpsi[100][1],psi_1.psi[1100][1],psi_1.fpGpsi[1100][1]);
	
		for(s_cntr=1;s_cntr<imx.s;++s_cntr)
		{

			for(j=0;j<s_cntr;++j)
			{
			   for(i=0;i<N;++i)
			   {
				
				//printf("chk %lf %lf %lf %lf\n",psi_1.psi[i][0],psi_1.fpGpsi[i][0],psi_1.psi[i][1],psi_1.fpGpsi[i][1]);
				 dbi = (double)(i);
				xv = xini + dx*dbi;
 				if(j==0)
				{
				


				
					c_psi_amp2 = psi_1.fpGpsi[i][0]*psi_1.fpGpsi[i][0] + psi_1.fpGpsi[i][1]*psi_1.fpGpsi[i][1]; 
				

					psi_1.ex_rhs(i,s_cntr-1);
					
					psi_1.im_rhs(i,s_cntr-1);
					

					psi_1.fpGpsi[i][0] = psi_1.psi[i][0] + dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[0][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[0][j][i];
					psi_1.fpGpsi[i][1] = psi_1.psi[i][1] + dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[1][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[1][j][i];

				

				}


				else
				{psi_1.fpGpsi[i][0]+=  dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[0][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[0][j][i];
				 psi_1.fpGpsi[i][1]+=  dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[1][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[1][j][i];

				 
				}



			   }


			}

		

			psi_1.do_forward_fft();
			


			for(i=0;i<N;++i)
			{
			
				
			 lambda = (k_grid[i]*k_grid[i])*imx.im_a[s_cntr][s_cntr]*dt;
			 kv = k_grid[i];		

			 psi_1.update_fft_fields(i,kv,lambda);
			 



			}

			


			psi_1.do_back_fft();
			




		}//RK stages concluded



		


		for(i=0;i<N;++i)
		{

			
			c_psi_amp2 = psi_1.fpGpsi[i][0]*psi_1.fpGpsi[i][0] + psi_1.fpGpsi[i][1]*psi_1.fpGpsi[i][1]; 
			
			psi_1.ex_rhs(i,imx.s-1);
			

			psi_1.im_rhs(i,imx.s-1);
			


			for(j=0;j<imx.s;++j)
			{	psi_1.psi[i][0]+=  dt*imx.ex_b[j]*psi_1.ex_K_psi[0][j][i]+ dt*imx.im_b[j]*psi_1.im_K_psi[0][j][i];
				psi_1.psi[i][1]+=  dt*imx.ex_b[j]*psi_1.ex_K_psi[1][j][i]+ dt*imx.im_b[j]*psi_1.im_K_psi[1][j][i];

				




			}

			psi_1.fpGpsi[i][0] = psi_1.psi[i][0];
			psi_1.fpGpsi[i][1] = psi_1.psi[i][1];

		


			if(isnan((psi_1.psi[i][0]+psi_1.psi[i][1])))
			{
				fail =1;
				//if(prt)
				printf("FAILED %d tcntr %d %lf %lf \n",i,tcntr,psi_1.psi[i][0],psi_1.psi[i][1]);

				break;

			}



		   
		}

		


	}///// ENd f Time Evolution /////////////////////////

	fclose(fp);
	fftw_cleanup_threads();
	*mass_err = psi_1.max_mass_err;
	


	if(prnt)
	printf("N %d\n Run mass los %lf\n",N,*(mass_err));


	
	return(psi_1.max_mass_err);








}






int main(int argc, char ** argv)
{

	double dt = 3e-4;
	double dx = 4e-3;
	double err[2],mass_loss;

	double dx_l=2e-3,dx_u = 4e-2;
	double dt_l= 1e-5,dt_u = 1e-2;

	double ddx = (dx_u-dx_l)/(20.0);
	double ddt = (dt_u-dt_l)/(20.0);

	FILE *fp = fopen("imex_ft.txt","w");
	

	//for(dt=dt_l;dt<=dt_u;dt+=ddt)
	{


		//for(dx = dx_l;dx<=dx_u;dx+=ddx)
		{
			dx = 2e-2;
			dt = 1e-3;
			mass_loss = run(dt,512,err,argc,argv,1,1);

			printf("%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),mass_loss,*err,*(err+1));
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),mass_loss,*err,*(err+1));
		}

	}

	fclose(fp);




}
