using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>

#define pie M_PI 



//////////GLobal constants/////////////

double m,n,T;

///////////////////////////////////////


//////////// ImEx RK Butcher Tableau /////
const int imex_s = 3;
const int im_s = imex_s;
const int ex_s = imex_s;
/*
double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  205.0/462.0,2.0/11.0,0.0,  2033.0/4620.0,21.0/110.0,2.0/11.0};
double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  5.0/6.0,0.0,0.0,  11.0/24.0,11.0/24.0,0.0};

double im_c[im_s] = {2.0/11.0,289.0/462.0,751.0/924.0};
double ex_c[ex_s] = {0.0,5.0/6.0,11.0/12.0};

double im_b[im_s] = {24.0/55.0,1.0/5.0,4.0/11.0};
double ex_b[ex_s] = {24.0/55.0,1.0/5.0,4.0/11.0};
*/

//	/	/	/	/	/	/	/	/	/	/	

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  41.0/154.0,2.0/11.0,0.0,  289.0/847.0,42.0/121.0,2.0/11.0};
double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  0.5,0.0,0.0,  0.5,0.5,0.0};

double im_c[im_s] = {2.0/11.0,69.0/154.0,67.0/77.0};
double ex_c[ex_s] = {0.0,0.5,1.0};

double im_b[im_s] = {1.0/3.0,1.0/3.0,1.0/3.0};
double ex_b[ex_s] = {1.0/3.0,1.0/3.0,1.0/3.0};


////////////////////////////////////////////


double b(double x,double Pv)
{
	double eps = 0.005;	

	if(x<=0.5)
	return(eps/(eps+Pv));
	else
	return(100.0*eps/(eps+Pv));


}

double f(double t,double x)
{
	double r;
	
	if(t==0.0)
	{
		
		r = ((double)rand())/((double) RAND_MAX);
		r = 0.8 + r*(1.2-0.8);
		//printf("%lf %lf\n",x,r);
	}

	

	else
	r=0.0;
	return(r);


}




double ex_vel(double t,double x,double Pv)
{
	double fval,bval,res;
	double rd = 1.0;
	
	fval = f(t,x);
	bval = b(x,Pv);

	return(fval+bval-rd*Pv);



}




void initialise(double *P,double *x,double *k,double dx,int N)
{

	int i,j; double di,dk;
	dk = 2.0*pie/(((double)(N))*dx);
	for(i=0;i<N;++i)
	{
	   di = (double)i;	
	   P[i] = 0.0;

	    if(i<=(N/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) N))*dk; 

		x[i] = 0.0+di*dx;

	//	printf("k %lf  %lf  dx %lf\n",dk,*(k+i),dx );

	}


}






int main()
{

	int N,t_steps;
	double box_len,dx,t_end,t_start,dt,diff;


///////////////////////////////File for data/////////////////////////////////////////////////////

	FILE *fp = fopen("data_pred.txt","w");
	FILE *fp2 = fopen("data2_pred.txt","w");
	
	FILE *fpmass = fopen("mass_pred.txt","w");

/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	m  = 1.0;
	n  = 1.0;
	T = 2.0*pie*m;
	diff = 0.02;
	N = 100;



/////////////////////////////////////////RK things/////////////////////////////////////////


	double im_K_P[ex_s][N],ex_K_P[im_s][N];


/////////////////////////////// Box & res. setting ///////////////////////////////////////////////

	box_len = 1.0;
	dx = box_len/(double(N-1));

////////////////////////////// Time & dt settings ////////////////////////////////////////////////
	
	t_start = 0.0;
	t_end = 10.0;
	t_steps = 20000;
	dt  = (t_end-t_start)/((double)t_steps);
	printf("dt %lf\n",dt);
	

////////////////////////////// Psi variables  /////////////////////////////////////////////////////

	double P[N],lap_val[2*N],lambda;

	
	double x[N],k_grid[N];
	
	fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_complex *K;
	fftw_complex *K_ft;

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	fftw_plan plan_imp_b;

	fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *N);

	K = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	K_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *N);


	plan_pois_f = fftw_plan_dft_1d(N, fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_pois_b = fftw_plan_dft_1d(N,fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);

	
	plan_imp_b = fftw_plan_dft_1d(N,K_ft, K, FFTW_BACKWARD, FFTW_ESTIMATE);

	initialise(P,x,k_grid,dx,N);

	int i,j;
	
	 	for(i=0;i<N;++i)
		{	
			
			fpGpsi[i][0] = P[i];
			fpGpsi[i][1] = 0.0 ;


		}


///////////////////////  Evolution ///////////////////////////////////////////////////////////////
	srand(time(0));

	int s_cntr,tcntr,printcntr,fail=0;
	printcntr = (int)(((double) t_steps)/100.0);	printf("%d\n",printcntr);
	if(t_steps<=100)
		printcntr = 1.0;
	double t,vel_val;



	

		
	

	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(5);t+=dt,++tcntr)
	{	
		
	
	
		  fftw_execute(plan_pois_f);

	 	for(i=0;i<N;++i)
		{	
			
			if((tcntr%printcntr)==0)  
			{
			  
			  
			  fprintf(fp,"%lf\t%lf\n",x[i],P[i]);
			   
			}
			
			lambda = k_grid[i]*k_grid[i]*im_a[0][0]*diff*dt;


			fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
			fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] + lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);

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


		for(s_cntr=1;s_cntr<imex_s;++s_cntr)
		{

			for(j=0;j<s_cntr;++j)
			{
			   for(i=0;i<N;++i)
			   {
				
				if(j==0)
				{

					
		    			vel_val = ex_vel(t+ex_c[s_cntr-1]*dt,x[i],fpGpsi[i][0]);

		    			ex_K_P[s_cntr-1][i] = vel_val; 
					im_K_P[s_cntr-1][i] = K[i][0]*diff; 

					fpGpsi[i][0] = P[i] + dt*ex_a[s_cntr][j]*ex_K_P[j][i]+ dt*im_a[s_cntr][j]*im_K_P[j][i];
					fpGpsi[i][1] = 0.0;

				}
			
				
				else
				{fpGpsi[i][0]+=  dt*ex_a[s_cntr][j]*ex_K_P[j][i]+ dt*im_a[s_cntr][j]*im_K_P[j][i];
				 fpGpsi[i][1]=  0.0;
				}
				
	
	
			   }

	
			}
	
			   fftw_execute(plan_pois_f);

			  for(i=0;i<N;++i)
			  {	
				lambda = k_grid[i]*k_grid[i]*im_a[s_cntr][s_cntr]*diff*dt;

				fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
				fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] + lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);

				

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
			  
			  // printf("%lf\n",t/t_end);
		}
			
		

		for(i=0;i<N;++i)
		{	
	
			


			vel_val = ex_vel(t+ex_c[imex_s-1]*dt,x[i],fpGpsi[i][0]);

		    	ex_K_P[imex_s-1][i] = vel_val;  
			im_K_P[imex_s-1][i] = K[i][0]*diff;  

	
							
			   for(j=0;j<imex_s;++j)
			   {	P[i]+=  dt*ex_b[j]*ex_K_P[j][i]+ dt*im_b[j]*im_K_P[j][i];
				
				
				if(isnan(P[i])||(P[i]<0.0))
			 	{

					if(isnan(P[i]))
					printf("!!!GONE NAN!!!!  %d\tt %lf\t%lf\n",tcntr,t,P[i]);
					if(P[i]<0.0)
					printf("!!!+tivity broken!!!  %d\tt %lf\t%lf\n",tcntr,t,P[i]);
					
					fail = 1;

					break;
				}
		
				
			    }

			if(fail)
			break;

			fpGpsi[i][0] = P[i];
			fpGpsi[i][1] = 0.0;


			


			

			if((tcntr%printcntr)==0)  
			fprintf(fp2,"%lf\t%lf\t%lf\n",x[i],P[i],fpGpsi[i][0]);

			

		   }



		if((tcntr%printcntr)==0)  
		{ fprintf(fp,"\n\n\n");
		  fprintf(fp2,"\n\n\n");
		  fprintf(fpmass,"%lf\n",t/t_end);
		  printf("dt %lf\tt %lf\n\n",dt,t);
		}
		


		  

	}///// ENd f Time Evolution /////////////////////////
	
	





	return(0);

}
