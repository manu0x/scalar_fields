

#include <stdio.h>
#include <math.h>



#include <fftw3.h>


#include <limits>


#include <stdlib.h>
#include <string.h>
#include <string>

#include <algorithm>
#include "../imex/imex_classes.cpp"
#include "../GPE/GPE_classes_3d_fid.cpp"

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


double run(double dt,int N,double *mass_err,int argc,char **argv,int prntfp,int prnt)
{

	int t_steps;
	double box_len,t,t_end,t_start,tk,dt_up,xval,dbi,dbj,dbk,dx;

	double x0[3];
	double *xv;
	double xvp[3];

	dt_up=dt;

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
	box_len = 1.0;
	x0[0]=-0.5*box_len;
	x0[1]=-0.5*box_len;
	x0[2]=-0.5*box_len;
	
	//n  = 2.0/box_len;
	N = 128;
	dx = box_len/(double(N));
	
	//N = ((int)(box_len/dx));
	int N3 = N*N*N;
	int N2 = N*N;
	printf("N2 %d\n",N2);

	xv = new double[3*N3];

////////////	Imex declare and data write file declare ////////////////////////



	char fp_name[30]("data_");


    int stages;
	char *imex_file; 
   


    stages = atoi(argv[1]);
	imex_file = argv[2];
    
	
	
	strcat(fp_name,imex_file);
	printf("ImEx table filename is %s and stages given by u is %d and out file %s \n",imex_file,stages,fp_name);

	FILE *fp = fopen(fp_name,"w");
	
	
	imex_table imx(stages);
    
	
  	imx.read_from_file(imex_file);
  	imx.print_table();
   
	

////////////////////////	Declare field /////////////////////////////////////////////
char *f1paramfile;

GPE_field_3d  psi_1(0,N,imx.s,16);

f1paramfile = argv[3];

psi_1.read_from_file(f1paramfile);
	

psi_1.print_params_set_kappa();







////////////////////////////// Time & da settings ////////////////////////////////////////////////

	t_start = psi_1.ti;
	t_end = 0.1*psi_1.t0;

	//dt = dt_up;*psi_1.aoft(t_start);
	
	t_steps = (int)((t_end-t_start)/dt);
	
	if(prnt)
	printf("dt %lf N %d\n",dt,N);
	printf("ti %lf t0 %lf %lf\n",psi_1.ti,psi_1.t0,psi_1.ti/psi_1.t0);


	
/////////////////////////////////////////RK things/////////////////////////////////////////
	int i,j;

////////////////////////////// k things /////////////

	double k_grid[N],ksqr;


//////////////////////////////  variables  /////////////////////////////////////////////////////



	double lambda;



    
	

	initialise_kgrid(k_grid,dx,N);

	//psi_1.initialise_random(k_grid) ; 
	//psi_1.read_from_initial();
	
	
	
	//printf("no of threads %d\n",omp_get_num_threads());
	//omp_set_num_threads(4);
	//printf("no of threads %d\n",omp_get_num_threads());
	


	


///////////////////////  Evolution ///////////////////////////////////////////////////////////////
    int ii,jj,kk;
	int s_cntr,tcntr,printcntr,fpcntr,err_cntr=1,fail=0;

	printcntr = (int)(((double) t_steps)/100.0);
	fpcntr = (int)(((double) t_steps)/10.0);
	if(t_steps<=100)
	{	printcntr = 1.0;
	
	}	//printf("%d\n",printcntr);
	double vel_val[2],sol[2],c_psi_amp2,delta,amp,avg_amp;
	
	double amp_ini,dsum,fac,lfac,vmax=-0.00000000000001;
	

	double kv[3];
	int ind[3];


	ii=-1;
	jj=-1;
	kk=0;

	psi_1.reset_consv_quant(1);

		for(i=0,ii=-1,jj=-1,kk=0;i<(N3);++i,++kk)
		{	

			

			if(i%N ==0)
            {
                kk=0;
                ++jj;
                if(i%N2==0)
                {
                    jj=0;
                    ++ii;


                }
			}
			

			dbi = (double)(ii); dbj = (double)(jj); dbk = (double)(kk);
			
			xvp[0] = x0[0] + dx*dbi; xvp[1] = x0[1] + dx*dbj; xvp[2] = x0[2] + dx*dbk;
			*(xv+3*i) = x0[0] + dx*dbi; *(xv+3*i+1) = x0[1] + dx*dbj; *(xv+3*i+2) = x0[2] + dx*dbk;
			psi_1.initialise_point(i,xvp);
	
		}


		if(prnt)
	printf("Initialization done\n");
	
	psi_1.set_field();
	
	
	printf("Starting Run..,\n");
////////////// START RUNNINGGGGGGGGGGG	///////////////////////////////////////////////////////////////////////
	
	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(30);t+=dt,++tcntr)
	{
		
		////////////////////////////	Adaptive step checks	/////////////////////////////
		
		
	/*	
		for ( i = 0; i < N3; i++)
		{
			if(vmax<fabs(psi_1.V_phi[i][0]))
			{			vmax = fabs(psi_1.V_phi[i][0]);
				
			}
		}
	

		if(vmax>(imx.ex_stb_r/dt))
		{		printf(" vmax dt %lf stb r %lf  %lf\n",vmax*dt,imx.ex_stb_r,fac);

			dt_up = 0.5*imx.ex_stb_r/vmax;
		}

		if(tcntr)
		dt = dt_up;*/

		//if(!tcntr)
		//dt = dt_up*psi_1.aoft(t);
		//////////////////////////////////////////////////////////////////////////////////////////
		if((tcntr%printcntr==0)&&prnt)
		printf("time %lf %lf vmax %lf  dt %lf\n",t/t_end,psi_1.mass,vmax,dt);
		psi_1.do_forward_fft();
		

		psi_1.reset_consv_quant();
		dsum = 0.0;
		
		
	 	for(i=0,ii=-1,jj=-1,kk=0;i<(N3);++i,++kk)
		{	

			

			if(i%N ==0)
            {
                kk=0;
                ++jj;
                if(i%N2==0)
                {
                    jj=0;
                    ++ii;


                }
			}
			

			dbi = (double)(ii); dbj = (double)(jj); dbk = (double)(kk);
			xvp[0] = x0[0] + dx*dbi; xvp[1] = x0[1] + dx*dbj; xvp[2] = x0[2] + dx*dbk;
			if(tcntr%err_cntr==0)
			{
				psi_1.cal_conserve_at_point(i,t,xvp, dx,!tcntr);
				

			}

			if((tcntr%fpcntr==0)&&(prntfp))
			{
				
				c_psi_amp2 = (psi_1.psi[i][0]*psi_1.psi[i][0]+psi_1.psi[i][1]*psi_1.psi[i][1]);
				delta = (c_psi_amp2)-1.0;

				dsum+=delta; 
				psi_1.sol(t,xvp,sol);
				fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",*(xv+3*i),*(xv+3*i+1),*(xv+3*i+2),psi_1.psi[i][0],psi_1.psi[i][1],sol[0],sol[1],c_psi_amp2,delta);


			}
			//tk = t+imx.im_c[0]*dt;
			//a = psi_1.aoft(tk);
			ksqr = k_grid[ii]*k_grid[ii] + k_grid[jj]*k_grid[jj] + k_grid[kk]*k_grid[kk];
			lambda = ksqr*imx.im_a[0][0]*psi_1.kppa*dt/(2.0);
			//printf("lmda  %lf\n",lambda);
			
			

			psi_1.update_fft_fields(i,ksqr,lambda);
			

		}
		if((tcntr%fpcntr==0)&&(prntfp))
			printf("t/t_end %lf dsum = %lf\n",t/t_end,dsum);
		if(tcntr%err_cntr==0)
		{	psi_1.conserve_err();
			

		}

		if((tcntr%fpcntr==0)&&(prntfp))
		fprintf(fp,"\n\n\n");
		
		
		psi_1.do_back_fft();
	
		

		for(s_cntr=1;s_cntr<imx.s;++s_cntr)
		{
			
			for(j=0;j<s_cntr;++j)
			{
			   for(i=0;i<N3;++i)
			   {
				
 				if(j==0)
				{
				
					tk = t+imx.ex_c[s_cntr-1]*dt;
					/*a = psi_1.aoft(tk);
					fac = a;*/

					

					lfac = psi_1.ex_rhs(i,s_cntr-1,tk,(xv+3*i));

					if(vmax<lfac)
					vmax = lfac;

				/*	if(psi_1.vmax*dt>imx.ex_stb_r)
					{
						printf("psi.vmax*dt  %lf stb_r %lf tcntr %d s_cntr-1 %d\n",psi_1.vmax*dt,imx.ex_stb_r,tcntr,s_cntr-1);
						//break;


					}

			*/
					//tk = t+imx.im_c[s_cntr-1]*dt;
					//a = psi_1.aoft(tk);
					psi_1.im_rhs(i,s_cntr-1);
		

					psi_1.fpGpsi[i][0] = psi_1.psi[i][0] + dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[0][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[0][j][i];
					

					

				}


				else
				{psi_1.fpGpsi[i][0]+=  dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[0][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[0][j][i];
				 psi_1.fpGpsi[i][1]+=  dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[1][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[1][j][i];
				}



			   }


			}

			if(vmax>(imx.ex_stb_r/dt))
			{printf("acntr %d  vmax da %lf stb r %lf\n",tcntr,vmax*dt,imx.ex_stb_r);
		//		da_up = 0.7*imx.ex_stb_r/vmax;
			}

			psi_1.do_forward_fft();
			
		


			for(i=0,ii=-1,jj=-1,kk=0;i<N3;++i,++kk)
			{
			
				if(i%N ==0)
            	{
                	kk=0;
                	++jj;
                	if(i%N2==0)
                	{
                    	jj=0;
                    	++ii;


                	}
				}
				
				
				//tk = t+imx.im_c[s_cntr]*dt;
				//a = psi_1.aoft(tk);
				ksqr = k_grid[ii]*k_grid[ii] + k_grid[jj]*k_grid[jj] + k_grid[kk]*k_grid[kk];
				
			 	lambda = ksqr*imx.im_a[s_cntr][s_cntr]*psi_1.kppa*dt/(2.0);
			 	

			 	psi_1.update_fft_fields(i,ksqr,lambda);
			 



			}

			


			psi_1.do_back_fft();
			



		}//RK stages concluded



		


		for(i=0,ii=-1,jj=-1,kk=0;i<N3;++i,++kk)
		{
			
				


			
			
			tk = t+imx.ex_c[imx.s-1]*dt;
			//a = psi_1.aoft(tk);
			psi_1.ex_rhs(i,imx.s-1,tk,(xv+3*i));

			
			//tk = t+imx.im_c[imx.s-1]*dt;
			//a = psi_1.aoft(tk);
			psi_1.im_rhs(i,imx.s-1);
		


			for(j=0;j<imx.s;++j)
			{	psi_1.psi[i][0]+=  dt*imx.ex_b[j]*psi_1.ex_K_psi[0][j][i]+ dt*imx.im_b[j]*psi_1.im_K_psi[0][j][i];
				psi_1.psi[i][1]+=  dt*imx.ex_b[j]*psi_1.ex_K_psi[1][j][i]+ dt*imx.im_b[j]*psi_1.im_K_psi[1][j][i];

				




			}

			psi_1.fpGpsi[i][0] = psi_1.psi[i][0];
			psi_1.fpGpsi[i][1] = psi_1.psi[i][1];

			


			if(isnan((psi_1.psi[i][0]+psi_1.psi[i][1])+(psi_1.V_phi[i][0])))
			{
				fail =1;
				//if(prt)
				printf("FAILED %d acntr %d %lf %lf %lf %lf\n",i,tcntr,psi_1.psi[i][0],psi_1.psi[i][1],psi_1.V_phi[i][0],psi_1.V_phi[i][1]);

				break;

			}



		   
		}

		


	}///// ENd f Time Evolution /////////////////////////

	fclose(fp);
	fftw_cleanup_threads();
	*mass_err = psi_1.max_mass_err;
	*(mass_err+1)= psi_1.max_mass_err;


	if(prnt)
	printf("N %d\n Run en los %lf  max sol err %lf\n",N,*(mass_err),psi_1.max_sol_err);


	
	return(psi_1.max_mass_err);








}






int main(int argc, char ** argv)
{

	double da = 3e-4;
	double dx = 4e-3;
	double mass_err[2],mass_loss;

	double dx_l=2e-3,dx_u = 4e-2;
	double da_l= 1e-5,da_u = 1e-2;

	double ddx = (dx_u-dx_l)/(20.0);
	double dda = (da_u-da_l)/(20.0);

	FILE *fp = fopen("imex_ft.txt","w");
	

	//for(da=da_l;da<=da_u;da+=dda)
	{


		//for(dx = dx_l;dx<=dx_u;dx+=ddx)
		{
			dx = 2e-2;
			da = 1.0e-4;
			mass_loss = run(da,512,mass_err,argc,argv,1,1);

			printf("%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,da,da/(dx*dx),mass_loss,*mass_err,*(mass_err+1));
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,da,da/(dx*dx),mass_loss,*mass_err,*(mass_err+1));
		}

	}

	fclose(fp);




}
