

#include <stdio.h>
#include <math.h>



#include <fftw3.h>


#include <limits>


#include <stdlib.h>
#include <string.h>
#include <string>

#include <algorithm>
#include "../imex/imex_classes.cpp"
#include "../GPE/GPE_classes_3d_exp_uni.cpp"

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


double run(double da,int N,double *mass_err,int argc,char **argv,int prntfp,int prnt)
{

	int a_steps;
	double box_len,a,a_end,a_start,ak,da_up,xval,dbi,dbj,dbk,dx;
	double HbyH0;
	double x0[3],xv[3];

	da_up=da;

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
	x0[0]=-0.5*box_len; x0[1]=-0.5*box_len; x0[2]=-0.5*box_len;
	
	//n  = 2.0/box_len;
	N = 128;
	dx = box_len/(double(N));
	
	//N = ((int)(box_len/dx));
	int N3 = N*N*N;
	int N2 = N*N;
	printf("N2 %d\n",N2);

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

	a_start = psi_1.ai;
	a_end = 1.0;
	
	a_steps = (int)((a_end-a_start)/da);
	
	if(prnt)
	printf("da %lf N %d\n",da,N);
	HbyH0= psi_1.HbyH0(a_start);
	printf("chking c %lf %lf\n",(a_start*a_start*HbyH0),psi_1.kppa/(a_start*a_start*a_start*HbyH0));


	
/////////////////////////////////////////RK things/////////////////////////////////////////
	int i,j;

////////////////////////////// k things /////////////

	double k_grid[N],ksqr;


//////////////////////////////  variables  /////////////////////////////////////////////////////



	double lambda;



    
	

	initialise_kgrid(k_grid,dx,N);

	//psi_1.initialise_random(k_grid) ; 
	psi_1.read_from_initial();
	
	
	
	//printf("no of threads %d\n",omp_get_num_threads());
	//omp_set_num_threads(4);
	//printf("no of threads %d\n",omp_get_num_threads());
	

	if(prnt)
	printf("Initialization done\n");
	
	psi_1.set_field();
	


///////////////////////  Evolution ///////////////////////////////////////////////////////////////
    int ii,jj,kk;
	int s_cntr,acntr,printcntr,fpcntr,err_cntr=1,fail=0;

	printcntr = (int)(((double) a_steps)/100.0);
	fpcntr = (int)(((double) a_steps)/10.0);
	if(a_steps<=100)
	{	printcntr = 1.0;
	
	}	//printf("%d\n",printcntr);
	double vel_val[2],c_psi[2],c_psi_amp2,delta,amp,avg_amp;
	
	double amp_ini,dsum,fac,vmax=-0.00000000000001;
	

	double kv[3];
	int ind[3];


	ii=-1;
	jj=-1;
	kk=0;

	psi_1.reset_consv_quant(1);
	

	
	
	printf("Starting Run..,\n");

	
	for(a=a_start,acntr=0;(a<=a_end)&&(!fail)&&(30);a+=da,++acntr)
	{
		
		////////////////////////////	Adaptive step checks	/////////////////////////////
		
		psi_1.solve_V(k_grid);

		HbyH0 = psi_1.HbyH0(a);
		fac = a*a*HbyH0;
		for ( i = 0; i < N3; i++)
		{
			if(vmax<fabs(psi_1.V_phi[i][0]/fac))
			{			vmax = fabs(psi_1.V_phi[i][0]/fac);
				
			}
		}
	

		if(vmax>(imx.ex_stb_r/da))
		{		printf(" vmax da %lf stb r %lf  %lf\n",vmax*da,imx.ex_stb_r,fac);

			//da_up = 0.9*imx.ex_stb_r/vmax;
		}



		da = da_up;
		//////////////////////////////////////////////////////////////////////////////////////////
		if((acntr%printcntr==0)&&prnt)
		printf("time %lf %lf %lf  da %lf\n",a/a_end,psi_1.mass,psi_1.mass,da);


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
			xv[0] = x0[0] + dx*dbi; xv[1] = x0[1] + dx*dbj; xv[2] = x0[2] + dx*dbk;
			if(acntr%err_cntr==0)
			{
				psi_1.cal_conserve_at_point(i, dx,!acntr);
				

			}

			if((acntr%fpcntr==0)&&(prntfp))
			{
				
				c_psi_amp2 = (psi_1.psi[i][0]*psi_1.psi[i][0]+psi_1.psi[i][1]*psi_1.psi[i][1]);
				delta = (c_psi_amp2/(3.0*psi_1.omega_m0) -1.0);
				dsum+=delta; 
				fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",xv[0],xv[1],xv[2],psi_1.psi[i][0],psi_1.psi[i][1],c_psi_amp2,delta);


			}
			ak = a+imx.im_c[0]*da;
			HbyH0 = psi_1.HbyH0(ak);
			ksqr = k_grid[ii]*k_grid[ii] + k_grid[jj]*k_grid[jj] + k_grid[kk]*k_grid[kk];
			lambda = ksqr*imx.im_a[0][0]*psi_1.kppa*da/(ak*ak*ak*HbyH0*2.0);
			//printf("lmda  %lf\n",lambda);
			
			

			psi_1.update_fft_fields(i,ksqr,lambda);
			

		}
		if((acntr%fpcntr==0)&&(prntfp))
			printf("a/a_end %lf dsum = %lf\n",a/a_end,dsum);
		if(acntr%err_cntr==0)
		{	psi_1.conserve_err();
			

		}

		if((acntr%fpcntr==0)&&(prntfp))
		fprintf(fp,"\n\n\n");
		
		
		psi_1.do_back_fft();
		psi_1.solve_V(k_grid);
		

		for(s_cntr=1;s_cntr<imx.s;++s_cntr)
		{
			
			for(j=0;j<s_cntr;++j)
			{
			   for(i=0;i<N3;++i)
			   {
				
 				if(j==0)
				{
				
					ak = a+imx.ex_c[s_cntr-1]*da;
					HbyH0 = psi_1.HbyH0(ak);
					fac = ak*ak*HbyH0;

					if(vmax<fabs(psi_1.V_phi[i][0]/fac))
					vmax = fabs(psi_1.V_phi[i][0]/fac);

					psi_1.ex_rhs(i,s_cntr-1,fac);
					
					
					
			
					ak = a+imx.im_c[s_cntr-1]*da;
					HbyH0 = psi_1.HbyH0(ak);
					psi_1.im_rhs(i,s_cntr-1,ak*ak*ak*HbyH0);

					
					

					
					

					psi_1.fpGpsi[i][0] = psi_1.psi[i][0] + da*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[0][j][i]+ da*imx.im_a[s_cntr][j]*psi_1.im_K_psi[0][j][i];
					psi_1.fpGpsi[i][1] = psi_1.psi[i][1] + da*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[1][j][i]+ da*imx.im_a[s_cntr][j]*psi_1.im_K_psi[1][j][i];
					

					

				}


				else
				{psi_1.fpGpsi[i][0]+=  da*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[0][j][i]+ da*imx.im_a[s_cntr][j]*psi_1.im_K_psi[0][j][i];
				 psi_1.fpGpsi[i][1]+=  da*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[1][j][i]+ da*imx.im_a[s_cntr][j]*psi_1.im_K_psi[1][j][i];
				}



			   }


			}

		//	if(vmax>(imx.ex_stb_r/da))
		//	{printf("acntr %d  vmax da %lf stb r %lf\n",acntr,vmax*da,imx.ex_stb_r);
		//		da_up = 0.7*imx.ex_stb_r/vmax;
		//	}

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
				
				
				ak = a+imx.im_c[s_cntr]*da;
				HbyH0 = psi_1.HbyH0(ak);
				ksqr = k_grid[ii]*k_grid[ii] + k_grid[jj]*k_grid[jj] + k_grid[kk]*k_grid[kk];
				
			 	lambda = ksqr*imx.im_a[s_cntr][s_cntr]*psi_1.kppa*da/(ak*ak*ak*HbyH0*2.0);
			 	

			 	psi_1.update_fft_fields(i,ksqr,lambda);
			 



			}

			


			psi_1.do_back_fft();
			psi_1.solve_V(k_grid);



		}//RK stages concluded



		


		for(i=0,ii=-1,jj=-1,kk=0;i<N3;++i,++kk)
		{
			
				


			
			
			ak = a+imx.ex_c[imx.s-1]*da;
			HbyH0 = psi_1.HbyH0(ak);
			psi_1.ex_rhs(i,imx.s-1,ak*ak*HbyH0);

			
			ak = a+imx.im_c[imx.s-1]*da;
			HbyH0 = psi_1.HbyH0(ak);
			psi_1.im_rhs(i,imx.s-1,ak*ak*ak*HbyH0);
		


			for(j=0;j<imx.s;++j)
			{	psi_1.psi[i][0]+=  da*imx.ex_b[j]*psi_1.ex_K_psi[0][j][i]+ da*imx.im_b[j]*psi_1.im_K_psi[0][j][i];
				psi_1.psi[i][1]+=  da*imx.ex_b[j]*psi_1.ex_K_psi[1][j][i]+ da*imx.im_b[j]*psi_1.im_K_psi[1][j][i];

				




			}

			psi_1.fpGpsi[i][0] = psi_1.psi[i][0];
			psi_1.fpGpsi[i][1] = psi_1.psi[i][1];

			


			if(isnan((psi_1.psi[i][0]+psi_1.psi[i][1])+(psi_1.V_phi[i][0])))
			{
				fail =1;
				//if(prt)
				printf("FAILED %d acntr %d %lf %lf %lf %lf\n",i,acntr,psi_1.psi[i][0],psi_1.psi[i][1],psi_1.V_phi[i][0],psi_1.V_phi[i][1]);

				break;

			}



		   
		}

		


	}///// ENd f Time Evolution /////////////////////////

	fclose(fp);
	fftw_cleanup_threads();
	*mass_err = psi_1.max_mass_err;
	*(mass_err+1)= psi_1.max_mass_err;


	if(prnt)
	printf("N %d\n Run en los %lf abs err %lf\n",N,*(mass_err),*(mass_err+1));


	
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
			da = 1e-4;
			mass_loss = run(da,512,mass_err,argc,argv,1,1);

			printf("%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,da,da/(dx*dx),mass_loss,*mass_err,*(mass_err+1));
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,da,da/(dx*dx),mass_loss,*mass_err,*(mass_err+1));
		}

	}

	fclose(fp);




}
