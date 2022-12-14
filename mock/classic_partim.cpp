using namespace std;

#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#define pie M_PI 

//////////GLobal constants/////////////

double m,n,T;

///////////////////////////////////////

double V(double psi_amp)
{

	return( 1.0 - 2.0*pie*pie*n/(m*m));


}


void vel(double v[2],double psi[2],double Vval,double lap_psi[2])
{

	v[0] = m*Vval*psi[1];
	v[1] = -m*Vval*psi[0];;



}


void initialise(fftw_complex *psi,double *k,double dx,int N)
{

	int i,j; double di,dk;
	dk = 1.0/(((double)(N-1))*dx);
	for(i=0;i<N;++i)
	{
	   di = (double)i;	
	   psi[i][0] = sin(2.0*pie*di*dx);
	   psi[i][1] = 0.0;

	    if(i<=(N/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) N))*dk; 

		//printf("k %lf  %lf\n",dk,*(k+i) );

	}


}






int main()
{

	int N,t_steps;
	double box_len,dx,t_end,t_start,dt;

///////////////////////////////////////////////////////////////////////////////////////////////

	double rk[3]={0.5,0.5,1.0};

///////////////////////////////File for data/////////////////////////////////////////////////////

	FILE *fp = fopen("data.txt","w");
	FILE *fpmass = fopen("mass.txt","w");

/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	m  = 1.0;
	n  = 10.0;
	T = 2.0*pie*m;
	N = 2000;

/////////////////////////////// Box & res. setting ///////////////////////////////////////////////

	box_len = 2.0;
	dx = box_len/(double(N-1));

////////////////////////////// Time & dt settings ////////////////////////////////////////////////
	
	t_start = 0.0;
	t_end = 2.0*T;
	t_steps = 100000;
	dt  = (t_end-t_start)/((double)t_steps);
	printf("dt %lf\n",dt);
	

////////////////////////////// Psi variables  /////////////////////////////////////////////////////

	double psi[N][2],lap_val[2*N],lambda;

	double k_grid[N];
	
	fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_complex *A;
	fftw_complex *A_ft;

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;

	fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *N);


	plan_pois_f = fftw_plan_dft_1d(N, fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_pois_b = fftw_plan_dft_1d(N,fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);

	initialise(fpGpsi,k_grid,dx,N);
///////////////////////  Evolution ///////////////////////////////////////////////////////////////

	int i,j,tcntr,printcntr,fail=0;
	printcntr = 100;//(int)(((double) t_steps)/100.0);	printf("%d\n",printcntr);
	double t,vel_val[2],c_psi[2],c_psi_amp,Vval,amp,avg_amp;



		
	

	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(tcntr<1000);t+=dt,++tcntr)
	{	

		avg_amp = 0.0;
		for(i=0;i<N;++i)
		{	



			if((tcntr%printcntr)==0  &&(!fail))  
			 {
				   
				   amp  = sqrt(fpGpsi[i][0]*fpGpsi[i][0] + fpGpsi[i][1]*fpGpsi[i][1]);
				   avg_amp+=amp;
				   fprintf(fp,"%lf\t%lf\t%lf\n",fpGpsi[i][0],fpGpsi[i][1],amp);
				   
			  }
				

				
				c_psi[0] = fpGpsi[i][0];
				c_psi[1] =  fpGpsi[i][1]; 
					
				psi[i][0] = fpGpsi[i][0];
				psi[i][1] =  fpGpsi[i][1]; 
				
				c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
				Vval = V(c_psi_amp);
				vel(vel_val,c_psi,Vval,&lap_val[2*i]);

				fpGpsi[i][0] = fpGpsi[i][0] + dt*vel_val[0];
				fpGpsi[i][1] = fpGpsi[i][1] + dt*vel_val[1];


				if(isnan(c_psi_amp))
				   {

					printf("%d failed at tcntr %d time %lf\n %lf %lf \n",i,tcntr,t/t_end,fpGpsi[i][0],fpGpsi[i][1]);
					fail =1;
					break;
					
				   }
	
	
		}

			
		fftw_execute(plan_pois_f);

		for(i=0;i<N;++i)
		{	
			lambda = 2.0*pie*2.0*pie*k_grid[i]*k_grid[i]*dt/(2.0*m);

			fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
			fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);

			fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)N);
			fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)N);

			//if(i==0)
			//fpGpsi_ft[i][0] = 0.0;
			//fpGpsi_ft[i][1] = 0.0;

		}

		fftw_execute(plan_pois_b);
		
///##################################################   Second Substep starts here    ####################################################################

		avg_amp = 0.0;
		for(i=0;i<N;++i)
		{	


				
				c_psi[0] = fpGpsi[i][0];
				c_psi[1] =  fpGpsi[i][1]; 
					
				
				
				c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
				Vval = V(c_psi_amp);
				vel(vel_val,c_psi,Vval,&lap_val[2*i]);

				fpGpsi[i][0] = psi[i][0] + dt*vel_val[0];
				fpGpsi[i][1] = psi[i][1] + dt*vel_val[1];
	
				   if(isnan(c_psi_amp))
				   {

					printf("%d failed at tcntr %d time %lf\n %lf %lf \n",i,tcntr,t/t_end,fpGpsi[i][0],fpGpsi[i][1]);
					fail =1;
					break;
					
				   }	
		}

			
		fftw_execute(plan_pois_f);

		for(i=0;i<N;++i)
		{	
			lambda = k_grid[i]*k_grid[i]*dt/(2.0*m);

			fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
			fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);

			fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)N);
			fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)N);

		}

		fftw_execute(plan_pois_b);
		

		if((tcntr%printcntr)==0)  
		{ fprintf(fp,"\n\n\n");
		  fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)(N)));
		}



			  

	}
	avg_amp = 0.0;
	for(i=0;i<N;++i)
	{
		amp  = sqrt(fpGpsi[i][0]*fpGpsi[i][0] + fpGpsi[i][1]*fpGpsi[i][1]);
		avg_amp+=amp;
		fprintf(fp,"%lf\t%lf\t%lf\n",fpGpsi[i][0],fpGpsi[i][1],amp);
		//fprintf(fp,"%lf\t%lf\t%lf\n",*(psi+2*i),*(psi+2*i+1),amp);

	}

	fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)N));

	





	return(0);

}
