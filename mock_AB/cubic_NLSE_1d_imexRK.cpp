

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

using namespace std;

#define pie M_PI



//////////GLobal constants/////////////

double bta,a,c,x0,T;

///////////////////////////////////////


//////////// ImEx RK Butcher Tableau /////
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

/*
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

double V(double psi_amp2)
{

	return( psi_amp2);


}




void ex_vel(double v[2],double psi[2],double Vval)
{

	v[0] = bta*Vval*psi[1];
	v[1] = -bta*Vval*psi[0];;



}


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


void initialise(fftw_complex *psi,double *k,double xini,double dx,int N,double ti=0.0)
{

	int i,j; double di,dk;
	double theta,x,fs;
	dk = 2.0*pie/(((double)(N))*dx);


	x = xini;
	for(i=0;i<N;++i,x+=dx)
	{
	   di = (double)i;
	   
	   theta = 0.5*c*(x-x0)-(0.25*c*c -a)*ti;
	   fs = (2.0*a/bta)/cosh(sqrt(a)*(x-x0-c*ti));
	   if(i<N)
	  { psi[i][0] = cos(theta)*fs;
	    psi[i][1] = sin(theta)*fs;
	  }
	   if(i<(N/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) N))*dk;



	}




}




double run(double dt,double dx,double *ens,double *stb_avg,int stb_any,int printfp,int prt)
{

	int N,t_steps;
	double box_len,t_end,t_start,xval,dbi,sol[2],stb_ini,xini,fs,x,theta;
	double sol_n2,err_n2,loc_ens;
	*ens = -1.0;;

///////////////////////////////File for data/////////////////////////////////////////////////////

	FILE *fp = fopen("data_ft.txt","w");
	FILE *fp2 = fopen("data2_ft.txt","w");
	FILE *fplap = fopen("lap_ft.txt","w");
	FILE *fpmass = fopen("mass_ft.txt","w");
	FILE *fptime = fopen("tm_ft.txt","w");

/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	bta  = -8.0;
	a = bta*bta/16.0;

	T = 5.0;







/////////////////////////////// Box & res. setting ///////////////////////////////////////////////

	box_len = 40.0;
	xini = -20.0;
	N=512;
	dx = box_len/(double(N));
	//N = ((int)(box_len/dx));

////////////////////////////// Time & dt settings ////////////////////////////////////////////////

	t_start = 0.0;
	t_end = T;
	if(stb_any)
	t_end = 0.01*t_end;
	t_steps = (int)((t_end-t_start)/dt);
	//dt  = (t_end-t_start)/((double)t_steps);
	if(prt)
	printf("dt %lf N %d\n",dt,N);



/////////////////////////////////////////RK things/////////////////////////////////////////


	double im_K_psi[2][ex_s][N],ex_K_psi[2][im_s][N];


////////////////////////////// Psi variables  /////////////////////////////////////////////////////




	double psi[N][2],lambda;

	


	double k_grid[N];

	fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_complex *K;
	fftw_complex *K_ft;

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	fftw_plan plan_imp_b;

	fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
	fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N));

	K = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
	K_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(N));


	plan_pois_f = fftw_plan_dft_1d(N, fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_pois_b = fftw_plan_dft_1d(N,fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);


	plan_imp_b = fftw_plan_dft_1d(N,K_ft, K, FFTW_BACKWARD, FFTW_ESTIMATE);

	initialise(psi,k_grid,xini,dx,N);
	

	int i,j;

	 	for(i=0;i<N;++i)
		{

			fpGpsi[i][0] = psi[i][0] ;
			fpGpsi[i][1] = psi[i][1] ;

	





		}


///////////////////////  Evolution ///////////////////////////////////////////////////////////////

	int s_cntr,tcntr,printcntr,fail=0;
	printcntr = (int)(((double) t_steps)/100.0);
	if(t_steps<=100)
		printcntr = 1.0;	//printf("%d\n",printcntr);
	double t,vel_val[2],c_psi[2],c_psi_amp2,Vval,amp2,avg_amp2;
	//double drc = 2.0*pie*n*2.0*pie*n;
	double fdt,amp2_ini;
	*stb_avg=0.0;





	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(1);t+=dt,++tcntr)
	{

		avg_amp2 = 0.0;

		  fftw_execute(plan_pois_f);

	 	for(i=0;i<(N);++i)
		{	



			if((tcntr%printcntr)==0||(stb_any))
			{
			  
			   amp2  = (psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1]);
			   avg_amp2+=amp2;
			   dbi = (double)(i);
			   x = xini+dbi*dx;
			   theta = 0.5*c*(x-x0)-(0.25*c*c -a)*t;
	   		   fs = (2.0*a/bta)/cosh(sqrt(a)*(x-x0-c*t));
			   sol[0] = cos(theta)*fs;
			   sol[1] = sin(theta)*fs;
			  if(printfp)
			  fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp2);
			  

			  if((i==(10))&&(printfp))
			   fprintf(fptime,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp2);
			}

			lambda = k_grid[i]*k_grid[i]*im_a[0][0]*dt;


			fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
			fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);

			fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)N);
			fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)N);

			K_ft[i][0] = -fpGpsi_ft[i][0]*(k_grid[i]*k_grid[i]);
			K_ft[i][1] = -fpGpsi_ft[i][1]*(k_grid[i]*k_grid[i]);

			if(i==0)
			{
			 K_ft[i][0] = 0.0;
			 K_ft[i][1] = 0.0;
			}

		}
		stb_ini = avg_amp2;
		if(t==t_start)
		amp2_ini = avg_amp2;

		fftw_execute(plan_pois_b);
		fftw_execute(plan_imp_b);



		for(s_cntr=1;s_cntr<imex_s;++s_cntr)
		{

			for(j=0;j<s_cntr;++j)
			{
			   for(i=0;i<N;++i)
			   {


				if(j==0)
				{
			
					

					c_psi[0] = fpGpsi[i][0];
		    		c_psi[1] =  fpGpsi[i][1];




		    		c_psi_amp2 = (c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
		    		Vval = V(c_psi_amp2);
		    		ex_vel(vel_val,c_psi,Vval);

		    		ex_K_psi[0][s_cntr-1][i] = vel_val[0];  ex_K_psi[1][s_cntr-1][i] = vel_val[1];
					im_K_psi[0][s_cntr-1][i] = -K[i][1];  im_K_psi[1][s_cntr-1][i] =  K[i][0];

					fpGpsi[i][0] = psi[i][0] + dt*ex_a[s_cntr][j]*ex_K_psi[0][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[0][j][i];
					fpGpsi[i][1] = psi[i][1] + dt*ex_a[s_cntr][j]*ex_K_psi[1][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[1][j][i];

				}


				else
				{fpGpsi[i][0]+=  dt*ex_a[s_cntr][j]*ex_K_psi[0][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[0][j][i];
				 fpGpsi[i][1]+=  dt*ex_a[s_cntr][j]*ex_K_psi[1][j][i]+ dt*im_a[s_cntr][j]*im_K_psi[1][j][i];
				}



			   }


			}

		

			   fftw_execute(plan_pois_f);

			  for(i=0;i<N;++i)
			  {
				lambda = k_grid[i]*k_grid[i]*im_a[s_cntr][s_cntr]*dt;

				fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
				fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);



				fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)N);
				fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)N);

				K_ft[i][0] = -fpGpsi_ft[i][0]*(k_grid[i]*k_grid[i]);
				K_ft[i][1] = -fpGpsi_ft[i][1]*(k_grid[i]*k_grid[i]);

				if(i==0)
				{
				 K_ft[i][0] = 0.0;
				 K_ft[i][1] = 0.0;
				}



			   }

			  fftw_execute(plan_pois_b);
			  fftw_execute(plan_imp_b);







		}//RK stages concluded

		if((tcntr%printcntr)==0)
		{
			  if(prt)
			   printf("%lf\t%lf\n",t/t_end,avg_amp2/((double)(N)));
		}

		avg_amp2=0.0;
		err_n2 = 0.0;
		sol_n2 = 0.0;


		for(i=0;i<N;++i)
		{

			c_psi[0] = fpGpsi[i][0];
		    	c_psi[1] =  fpGpsi[i][1];

			c_psi_amp2 = (c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
		    	Vval = V(c_psi_amp2);
		    	ex_vel(vel_val,c_psi,Vval);

		    	ex_K_psi[0][imex_s-1][i] = vel_val[0];  ex_K_psi[1][imex_s-1][i] = vel_val[1];
			im_K_psi[0][imex_s-1][i] = -K[i][1];  im_K_psi[1][imex_s-1][i] =  K[i][0];


			   for(j=0;j<imex_s;++j)
			   {	psi[i][0]+=  dt*ex_b[j]*ex_K_psi[0][j][i]+ dt*im_b[j]*im_K_psi[0][j][i];
				    psi[i][1]+=  dt*ex_b[j]*ex_K_psi[1][j][i]+ dt*im_b[j]*im_K_psi[1][j][i];

			   }


			fpGpsi[i][0] = psi[i][0];
			fpGpsi[i][1] = psi[i][1];


			if(isnan((psi[i][0])+psi[i][1]))
			{
				fail =1;
				if(prt)
				printf("FAILED %d tcntr %d %lf %lf\n",i,tcntr,psi[i][0],psi[i][0]);

				break;

			}


			dbi = (double)(i);
			x = xini+dbi*dx;
			theta = 0.5*c*(x-x0)-(0.25*c*c -a)*(t+dt);
	   		fs = (2.0*a/bta)/cosh(sqrt(a)*(x-x0-c*(t+dt)));
			sol[0] = cos(theta)*fs;
			sol[1] = sin(theta)*fs;

			(err_n2)+=((sol[0]-psi[i][0])*(sol[0]-psi[i][0])+(sol[1]-psi[i][1])*(sol[1]-psi[i][1]));
			sol_n2+= (sol[0]*sol[0] + sol[1]*sol[1]);





			amp2  = (psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1]);
			avg_amp2+=amp2;

			if((tcntr%printcntr)==0)
			{
			  if(printfp)
			  fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp2);
			}


		}

		err_n2 = sqrt(err_n2);
		sol_n2 = sqrt(sol_n2);

		loc_ens = err_n2/sol_n2;
		if(loc_ens>(*ens))
			*ens = loc_ens;


		(*stb_avg)+=(avg_amp2/stb_ini);


		



		if((stb_any)&&(fail))
		{

			err_n2 = 1e5;
			*stb_avg = std::numeric_limits<double>::infinity();
			return(1e5);

		}
		else
		if(((100.0*fabs(avg_amp2-amp2_ini)/amp2_ini)>=1e3)||(fail)||((loc_ens)>=1e3))
		{

			err_n2 = -1e3;
			if((100.0*fabs(avg_amp2-amp2_ini)/amp2_ini)>=1e3)
			printf("Mass blow up %lf\n",(100.0*fabs(avg_amp2-amp2_ini)/amp2_ini));
			if((loc_ens)>=1e3)
			printf("Absolute err blow up %lf\n",loc_ens);
			return(-1e3);

		}



		if((tcntr%printcntr)==0)
		{
		  if(printfp)
		  {fprintf(fp,"\n\n\n");
		   fprintf(fp2,"\n\n\n");
		   fprintf(fpmass,"%lf\t%lf\t%lf\n",t/t_end,fabs(avg_amp2-amp2_ini)/amp2_ini,loc_ens);
		  }
		  if(prt)
		  printf("%lf\t%.10lf\n\n",t/t_end,avg_amp2/((double)(N)));
		}





	}///// ENd f Time Evolution /////////////////////////
	*stb_avg = (*stb_avg)/((double) tcntr);
	avg_amp2 = 0.0;
	for(i=0;i<N;++i)
	{	amp2  = (fpGpsi[i][0]*fpGpsi[i][0] + fpGpsi[i][1]*fpGpsi[i][1]);
		avg_amp2+=amp2;

	}
	if(printfp)
	fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp2/((double)N));


	if(prt)
	printf("N %d\n Mass loss %lf Ens %lf\n",N,fabs(avg_amp2-amp2_ini)/amp2_ini,*ens);



	return(fabs(avg_amp2-amp2_ini)/amp2_ini);









}






int main()
{

	double dt = 3e-4;
	double dx = 4e-3;
	double ens,m_loss,stb_avg;

	double dx_l=1e-3,dx_u = 4e-2;
	double dt_l= 1e-5 ,dt_u = 1e-1;

	double ddx = (dx_u-dx_l)/(6.0);
	double ddt = (dt_u-dt_l)/(6.0);

	FILE *fp = fopen("imex_ft.txt","w");

imex_table ch(3);
 char nn[10] = "test";
 
  ch.read_from_file(nn);
  ch.print_table();

	for(dt=dt_l;dt<=dt_u;dt+=ddt)
	{


		//for(dx = dx_l;dx<=dx_u;dx+=ddx)
		{
			//dx = 1e-3;
			//dt = 1e-4;
			//m_loss = run(dt,dx,&ens,&stb_avg,0,0,0);

			//printf("%lf\t%lf\t%lf\t%.10lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),ens,m_loss,stb_avg);
			//fprintf(fp,"%lf\t%lf\t%lf\t%.10lf\t%.10lf\t%.10lf\n",dx,dt,dt/(dx*dx),ens,m_loss,stb_avg);
		}

	}





}
