

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


/*//////////// ImEx RK Butcher Tableau /////
const int imex_s = 3;
const int im_s = imex_s;
const int ex_s = imex_s;

//	/	/	/	/ SCHEME 17	/	/	/	/	/	/

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  205.0/462.0,2.0/11.0,0.0,  2033.0/4620.0,21.0/110.0,2.0/11.0};
double im_c[im_s] = {2.0/11.0,289.0/462.0,751.0/924.0};
double im_b[im_s] = {24.0/55.0,1.0/5.0,4.0/11.0};

double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  5.0/6.0,0.0,0.0,  11.0/24.0,11.0/24.0,0.0};
double ex_c[ex_s] = {0.0,5.0/6.0,11.0/12.0};
double ex_b[ex_s] = {24.0/55.0,1.0/5.0,4.0/11.0};


//	/	/	/	/ SCHEME 20	/	/	/	/	/	/

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  41.0/154.0,2.0/11.0,0.0,  289.0/847.0,42.0/121.0,2.0/11.0};
double im_c[im_s] = {2.0/11.0,69.0/154.0,67.0/77.0};
double im_b[im_s] = {1.0/3.0,1.0/3.0,1.0/3.0};

double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  0.5,0.0,0.0,  0.5,0.5,0.0};
double ex_c[ex_s] = {0.0,0.5,1.0};
double ex_b[ex_s] = {1.0/3.0,1.0/3.0,1.0/3.0};


////////////////////////////////////////////


//	/	/	/	/ SCHEME 22	/	/	/	/	/	/

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  2829.0/9317.0,2.0/11.0,0.0,  148529.0/428582.0,7.0/23.0,2.0/11.0};
double im_c[im_s] = {2.0/11.0,4523.0/9317.0,15517.0/18634.0};
double im_b[im_s] = {1.0/3.0,1.0/3.0,1.0/3.0};

double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  0.5,0.0,0.0,  0.5,0.5,0.0};
double ex_c[ex_s] = {0.0,0.5,1.0};
double ex_b[ex_s] = {1.0/3.0,1.0/3.0,1.0/3.0};


////////////////////////////////////////////

//	/	/	/	/ SCHEME 23	/	/	/	/	/	/

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  2583.0/13310.0,2.0/11.0,0.0,  39731.0/139755.0,10.0/21.0,2.0/11.0};
double im_c[im_s] = {2.0/11.0,5003.0/13310.0,6271.0/6655.0};
double im_b[im_s] = {1.0/3.0,1.0/3.0,1.0/3.0};

double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  0.5,0.0,0.0,  0.5,0.5,0.0};
double ex_c[ex_s] = {0.0,0.5,1.0};
double ex_b[ex_s] = {1.0/3.0,1.0/3.0,1.0/3.0};

*/
////////////////////////////////////////////



void lap(double *lpsi,fftw_complex *psi, double dx,int N)
{
	int l1,l2,r1,r2,i;

	double vl1,vl2,vr1,vr2,vc;
   for(i=0;i<(N-1);++i)
   {	l1 = ((N-1)+ (i-1))%(N-1);
	l2 = ((N-1)+ (i-2))%(N-1);

	r1 = (i+1)%(N-1);
	r2 = (i+2)%(N-1);


	vl1 = psi[l1][0];
	vl2 = psi[l2][0];

	vr1 = psi[r1][0];
	vr2 = psi[r2][0];

	vc = psi[i][0];

	*(lpsi+2*i) =  (-vr2 + 16.0*vr1 - 30.0*vc + 16.0*vl1 - vl2)/(12.0*dx*dx);

	//if(i==0||i==1||i==2||i==(N-1))
	//{printf("\n%d\t%d\t%d\t%d\t%d\n",l2,l1,i,r1,r2);
	// printf("%lf\t%lf\t%lf\t%lf\t%lf\n",vl2,vl1,vc,vr1,vr2);
	// printf("%lf\t%lf\n\n",(-vr2 + 16.0*vr1 - 30.0*vc + 16.0*vl1 - vl2),*(lpsi+2*i));
	//}


	vl1 = psi[l1][1];
	vl2 = psi[l2][1];

	vr1 = psi[r1][1];
	vr2 = psi[r2][1];

	vc = psi[i][1];
	*(lpsi+2*i+1) = (-vr2 + 16.0*vr1 - 30.0*vc + 16.0*vl1 - vl2)/(12.0*dx*dx);

    }

	*(lpsi+2*i) = *(lpsi);
	*(lpsi+2*i+1) = *(lpsi+1);

}


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




double run(double dt,double dx,double *abs_err,int argc,char **argv,double *stb_avg,int stb_any,int printfp,int prt)
{

	int N,t_steps;
	double box_len,t_end,t_start,xval,dbi,dbj,dbk;
	double x0[2],xv[2];

	
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
	x0[0]=-10.0; x0[1]=10.0;
	box_len = x0[1]-x0[0];
	//n  = 2.0/box_len;
	N = 512;
	dx = box_len/(double(N));
	
	//N = ((int)(box_len/dx));
	//int N3 = N*N*N;
	int N2 = N*N;
	printf("N2 %d\n",N2);
////////////////////////////// Time & dt settings ////////////////////////////////////////////////

	t_start = 0.0;
	t_end = 10.0;
	
	t_steps = (int)((t_end-t_start)/dt);
	
	if(prt)
	printf("dt %lf N %d\n",dt,N);


	
/////////////////////////////////////////RK things/////////////////////////////////////////
	int i,j;




////////////////////////////// Psi variables  /////////////////////////////////////////////////////



	double lambda;


	char fp_name[30]("imex_ft_");

	char fp_phi1_r[20]("phi1_r_ini.txt");
	char fp_phi1_i[20]("phi1_i_ini.txt");

	char fp_phi2_r[20]("phi2_r_ini.txt");
	char fp_phi2_i[20]("phi2_i_ini.txt");
	
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


    GPE_field_2d  psi_1(0,N,3,8);
    GPE_field_2d  psi_2(1,N,3,8);

    psi_1.read_from_file(f1paramfile);
	psi_2.read_from_file(f2paramfile);

    psi_1.print_params();
	psi_2.print_params();

	psi_1.initialise_from_file(fp_phi1_r,fp_phi1_i); 
	psi_2.initialise_from_file(fp_phi2_r,fp_phi2_i); 
	//psi_2.initialise_from_file();

	double k_grid[N];

	initialise(k_grid,dx,N);

	
	psi_1.set_field();
	psi_2.set_field();


///////////////////////  Evolution ///////////////////////////////////////////////////////////////
    int ii,jj,kk;
	int s_cntr,tcntr,printcntr,fpcntr,fail=0;

	printcntr = (int)(((double) t_steps)/100.0);
	fpcntr = (int)(((double) t_steps)/10.0);
	if(t_steps<=100)
		printcntr = 1.0;	//printf("%d\n",printcntr);
	double t,vel_val[2],c_psi[2],c_psi_amp2[2],Vval,amp,avg_amp;
	
	double fdt,amp_ini;
	*stb_avg=0.0;

	double kv[2],Rv[2][2];


	ii=-1;
	jj=-1;
	kk=0;


	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&( 1);t+=dt,++tcntr)
	{

		avg_amp = 0.0;

		psi_1.do_forward_fft();
		psi_2.do_forward_fft();

	 	for(i=0,ii=-1,jj=0;i<(N2);++i,++jj)
		{	

			if((i%N)==0)
			{   jj=0;
				++ii;
				
			}
		//	printf("i,j,k  %d %d %d\n",ii,jj,kk);
	
			lambda = (k_grid[ii]*k_grid[ii]+k_grid[jj]*k_grid[jj])*imx.im_a[0][0]*dt/(2.0);
			kv[0] = k_grid[ii];		kv[1] = k_grid[jj];

			psi_1.update_fft_fields(i,kv,lambda);
			psi_2.update_fft_fields(i,kv,lambda);

		}
		

		psi_1.do_back_fft();
		psi_2.do_back_fft();

		/*
		for(i=0;i<(N);++i)
		{	dbi = (double)i;
			sol[0] = cos(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			sol[1] = -sin(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);


			fdt = (fpGpsi[i-1][0]+fpGpsi[i+1][0]-2.0*fpGpsi[i][0])/(dx*dx);
			if(i==0)
			fdt = (fpGpsi[N-1][0]+fpGpsi[i+1][0]-2.0*fpGpsi[i][0])/(dx*dx);
			if(i==(N-1))
			fdt = (fpGpsi[i-1][0]+fpGpsi[0][0]-2.0*fpGpsi[i][0])/(dx*dx);
			fprintf(fplap,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,K[i][0],K[i][1],-drc*fpGpsi[i][0],-drc*fpGpsi[i][1],fdt);

		}

		//printf("l do\n");
		*/

		for(s_cntr=1;s_cntr<imx.s;++s_cntr)
		{

			for(j=0;j<s_cntr;++j)
			{
			   for(i=0,ii=-1,jj=0;i<N2;++i,++jj)
			   {
				if((i%N)==0)
			    {   jj=0;
					++ii;
			
				}

				 dbi = (double)(ii); dbj = (double)(jj);
				xv[0] = x0[0] + dx*dbi; xv[1] = x0[1] + dx*dbj;
 				if(j==0)
				{
				/*	if(s_cntr==1)
					{
					  c_psi[0] = psi[i][0];
		    			  c_psi[1] =  psi[i][1];


					}
					else
					*/


					Rv[0][0] = psi_2.R[i][0];  Rv[0][1] = psi_2.R[i][1];
					Rv[1][0] = psi_1.R[i][0];  Rv[1][1] = psi_1.R[i][1];
					c_psi_amp2[0] = psi_1.fpGpsi[i][0]*psi_1.fpGpsi[i][0] + psi_1.fpGpsi[i][1]*psi_1.fpGpsi[i][1]; 
					c_psi_amp2[1] = psi_2.fpGpsi[i][0]*psi_2.fpGpsi[i][0] + psi_2.fpGpsi[i][1]*psi_2.fpGpsi[i][1]; 

					psi_1.ex_rhs(i,s_cntr-1,c_psi_amp2[0],c_psi_amp2[1],Rv[0],xv);
					psi_2.ex_rhs(i,s_cntr-1,c_psi_amp2[1],c_psi_amp2[0],Rv[1],xv);

					psi_1.im_rhs(i,s_cntr-1);
					psi_2.im_rhs(i,s_cntr-1);

					psi_1.fpGpsi[i][0] = psi_1.psi[i][0] + dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[0][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[0][j][i];
					psi_1.fpGpsi[i][1] = psi_1.psi[i][1] + dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[1][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[1][j][i];

					psi_2.fpGpsi[i][0] = psi_1.psi[i][0] + dt*imx.ex_a[s_cntr][j]*psi_2.ex_K_psi[0][j][i]+ dt*imx.im_a[s_cntr][j]*psi_2.im_K_psi[0][j][i];
					psi_2.fpGpsi[i][1] = psi_1.psi[i][1] + dt*imx.ex_a[s_cntr][j]*psi_2.ex_K_psi[1][j][i]+ dt*imx.im_a[s_cntr][j]*psi_2.im_K_psi[1][j][i];

				}


				else
				{psi_1.fpGpsi[i][0]+=  dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[0][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[0][j][i];
				 psi_1.fpGpsi[i][1]+=  dt*imx.ex_a[s_cntr][j]*psi_1.ex_K_psi[1][j][i]+ dt*imx.im_a[s_cntr][j]*psi_1.im_K_psi[1][j][i];

				 psi_2.fpGpsi[i][0]+=  dt*imx.ex_a[s_cntr][j]*psi_2.ex_K_psi[0][j][i]+ dt*imx.im_a[s_cntr][j]*psi_2.im_K_psi[0][j][i];
				 psi_2.fpGpsi[i][1]+=  dt*imx.ex_a[s_cntr][j]*psi_2.ex_K_psi[1][j][i]+ dt*imx.im_a[s_cntr][j]*psi_2.im_K_psi[1][j][i];
				}



			   }


			}

		

			psi_1.do_forward_fft();
			psi_2.do_forward_fft();


			for(i=0,ii=-1,jj=0;i<N2;++i,++jj)
			{
			
				if((i%N)==0)
			    {   jj=0;
					++ii;
			
				}
			 lambda = (k_grid[ii]*k_grid[ii]+k_grid[jj]*k_grid[jj])*imx.im_a[s_cntr][s_cntr]*dt/(2.0);
			 kv[0] = k_grid[ii];		kv[1] = k_grid[jj];

			 psi_1.update_fft_fields(i,kv,lambda);
			 psi_2.update_fft_fields(i,kv,lambda);



			}

			


			psi_1.do_back_fft();
			psi_2.do_back_fft();




		}//RK stages concluded

		if((tcntr%printcntr)==0)
		{
			  if(prt)
			   printf("%lf\t%lf\n",t/t_end,avg_amp/((double)(N2)));
		}

		avg_amp=0.0;
		*abs_err = 0.0;


		for(i=0,ii=-1,jj=0;i<N2;++i,++jj)
		{

			if((i%N)==0)
			{   jj=0;
				++ii;
			
			}

			Rv[0][0] = psi_2.R[i][0];  Rv[0][1] = psi_2.R[i][1];
			Rv[1][0] = psi_1.R[i][0];  Rv[1][1] = psi_1.R[i][1];
			c_psi_amp2[0] = psi_1.fpGpsi[i][0]*psi_1.fpGpsi[i][0] + psi_1.fpGpsi[i][1]*psi_1.fpGpsi[i][1]; 
			c_psi_amp2[1] = psi_2.fpGpsi[i][0]*psi_2.fpGpsi[i][0] + psi_2.fpGpsi[i][1]*psi_2.fpGpsi[i][1]; 

			psi_1.ex_rhs(i,imx.s-1,c_psi_amp2[0],c_psi_amp2[1],Rv[0],xv);
			psi_2.ex_rhs(i,imx.s-1,c_psi_amp2[1],c_psi_amp2[0],Rv[1],xv);

			psi_1.im_rhs(i,imx.s-1);
			psi_2.im_rhs(i,imx.s-1);


			for(j=0;j<imx.s;++j)
			{	psi_1.psi[i][0]+=  dt*imx.ex_b[j]*psi_1.ex_K_psi[0][j][i]+ dt*imx.im_b[j]*psi_1.im_K_psi[0][j][i];
				psi_1.psi[i][1]+=  dt*imx.ex_b[j]*psi_1.ex_K_psi[1][j][i]+ dt*imx.im_b[j]*psi_1.im_K_psi[1][j][i];

				psi_2.psi[i][0]+=  dt*imx.ex_b[j]*psi_2.ex_K_psi[0][j][i]+ dt*imx.im_b[j]*psi_2.im_K_psi[0][j][i];
				psi_2.psi[i][1]+=  dt*imx.ex_b[j]*psi_2.ex_K_psi[1][j][i]+ dt*imx.im_b[j]*psi_2.im_K_psi[1][j][i];




			}

			psi_1.fpGpsi[i][0] = psi_1.psi[i][0];
			psi_1.fpGpsi[i][1] = psi_1.psi[i][1];

			psi_2.fpGpsi[i][0] = psi_2.psi[i][0];
			psi_2.fpGpsi[i][1] = psi_2.psi[i][1];


			if(isnan((psi_1.psi[i][0]+psi_1.psi[i][1])+(psi_2.psi[i][0]+psi_2.psi[i][1])))
			{
				fail =1;
				if(prt)
				printf("FAILED %d tcntr %d %lf %lf\n",i,tcntr,psi_1.psi[i][0],psi_1.psi[i][0]);

				break;

			}


	






			

			if((tcntr%printcntr)==0)
			{
			 // if(printfp)
			 // fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp);
			}


		   }

			














	}///// ENd f Time Evolution /////////////////////////

	//if(printfp)
	//fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)N2));


	if(prt)
	printf("N %d\n Run en los %lf abs err %lf\n",N,100.0*fabs(avg_amp-amp_ini)/amp_ini,*abs_err);



	return(0.0);








}






int main(int argc, char ** argv)
{

	double dt = 3e-4;
	double dx = 4e-3;
	double abs_err,en_loss,stb_avg;

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
			dt = 1e-4;
			en_loss = run(dt,dx,&abs_err,argc,argv,&stb_avg,0,1,1);

			printf("%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err,stb_avg);
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err,stb_avg);
		}

	}





}
