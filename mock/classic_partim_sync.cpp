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



void vel_all(double v[2],double psi[2],double Vval,double lap_psi[2])
{

	v[0] = 2.0*pie*(Vval*psi[1] - lap_psi[1]/(2.0*m*m))/T;
	v[1] = 2.0*pie*(-Vval*psi[0] + lap_psi[0]/(2.0*m*m))/T;



}

void vel(double v[2],double psi[2],double Vval,double lap_psi[2])
{

	v[0] = m*Vval*psi[1];
	v[1] = -m*Vval*psi[0];;



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


void initialise(fftw_complex *psi,double *k,double dx,int N)
{

	int i,j; double di,dk;
	dk = 1.0/(((double)(N-1))*dx);
	for(i=0;i<N;++i)
	{
	   di = (double)i;	
	   psi[i][0] = sin(2.0*pie*n*di*dx);
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
	FILE *fplap = fopen("lap.txt","w");
	FILE *fpmass = fopen("mass.txt","w");

/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	m  = 1.0;
	n  = 1.0;
	T = 2.0*pie*m;
	N = 2000;

/////////////////////////////// Box & res. setting ///////////////////////////////////////////////

	box_len = 2.0;
	dx = box_len/(double(N-1));

////////////////////////////// Time & dt settings ////////////////////////////////////////////////
	
	t_start = 0.0;
	t_end = 2.0*T;
	t_steps = 200000;
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
	printcntr = (int)(((double) t_steps)/100.0);	printf("%d\n",printcntr);
	double t,vel_val[2],c_psi[2],c_psi_amp,Vval,amp,avg_amp;



		
	

	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(1);t+=dt,++tcntr)
	{	
		lap(lap_val,fpGpsi,dx,N);
		avg_amp = 0.0;
		for(i=0;i<N;++i)
		{	



			if((tcntr%printcntr)==0  &&(!fail))  
			 {
				   
				   amp  = sqrt(fpGpsi[i][0]*fpGpsi[i][0] + fpGpsi[i][1]*fpGpsi[i][1]);
				   avg_amp+=amp;
				   fprintf(fp,"%lf\t%lf\t%lf\n",fpGpsi[i][0],fpGpsi[i][1],amp);
				   fprintf(fplap,"%lf\t%lf\t%lf\t%lf\n",i*dx,lap_val[2*i],lap_val[2*i+1],-(2.0*pie)*(2.0*pie)*sin(2.0*pie*i*dx));
				   
			  }
				

				
				c_psi[0] = fpGpsi[i][0];
				c_psi[1] =  fpGpsi[i][1]; 
					
				psi[i][0] = fpGpsi[i][0];
				psi[i][1] =  fpGpsi[i][1]; 
				
				c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
				Vval = V(c_psi_amp);
				vel_all(vel_val,c_psi,Vval,&lap_val[2*i]);

				fpGpsi[i][0] = fpGpsi[i][0] ;//+ dt*vel_val[0];
				fpGpsi[i][1] = fpGpsi[i][1] ;//+ dt*vel_val[1];


				if(isnan(c_psi_amp))
				   {

					printf("%d failed at tcntr %d time %lf\n %lf %lf \n",i,tcntr,t/t_end,fpGpsi[i][0],fpGpsi[i][1]);
					fail =1;
					break;
					
				   }
	
	
		}

			

		
///##################################################   Second Substep starts here    ####################################################################

	
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
