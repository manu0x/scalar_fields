using namespace std;

#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#define pie M_PI 



//////////GLobal constants/////////////

double m,n,T;

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

double V(double psi_amp)
{

	return( 1.0 - 2.0*pie*pie*n*n/(m*m));


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

void ex_vel(double v[2],double psi[2],double Vval)
{

	v[0] = m*Vval*psi[1];
	v[1] = -m*Vval*psi[0];;



}

void simple_lap(double *lpsi,fftw_complex *psi, double dx,int N)
{
	int l1,l2,r1,r2,i;

	double vl1,vl2,vr1,vr2,vc;
   for(i=1;i<(N-1);++i)
   {	l1 = i-1;
	r1 = (i+1);

	if(i==0)
	l1 = N-2;
	if(i==0)
	r1 = N+2;
	vl1 = psi[l1][0];
	

	vr1 = psi[r1][0];
	

	vc = psi[i][0];
	
	*(lpsi+2*i) =  (vl1+vr1-2.0*vc)/(dx*dx);

	if(i==0)
	
	printf("lap i 0 %lf %lf %lf %d %d\n\n",*(lpsi+2*i),vl1,vr1,l1,r1);
	//if(i==0||i==1||i==2||i==(N-1))
	//{printf("\n%d\t%d\t%d\t%d\t%d\n",l2,l1,i,r1,r2);
	// printf("%lf\t%lf\t%lf\t%lf\t%lf\n",vl2,vl1,vc,vr1,vr2);
	// printf("%lf\t%lf\n\n",(-vr2 + 16.0*vr1 - 30.0*vc + 16.0*vl1 - vl2),*(lpsi+2*i));
	//}


	vl1 = psi[l1][1];
	

	vr1 = psi[r1][1];
	

	vc = psi[i][1];
	
	*(lpsi+2*i+1) =  (vl1+vr1-2.0*vc)/(dx*dx);

    }		
	
	
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


void initialise(fftw_complex *psi,double *k,double dx,int N,int spd)
{

	int i,j; double di,dk;
	dk = 2.0*pie/(((double)(N))*dx);
	int nN = N+2*spd;
	for(i=0;i<N;++i)
	{
	   di = (double)i;	
	   psi[i][0] = sin(2.0*pie*n*di*dx);
	   psi[i][1] = 0.0;

           
	  
	    if(i<=(N/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) N))*dk; 

	//	printf("k %lf  %lf  dx %lf\n",dk,*(k+i),dx );


	}

	


}






int main()
{

	int N,t_steps;
	double box_len,dx,t_end,t_start,dt,xval,dbi,sol[2];


///////////////////////////////File for data/////////////////////////////////////////////////////

	FILE *fp = fopen("data.txt","w");
	FILE *fp2 = fopen("data2.txt","w");
	FILE *fplap = fopen("lap.txt","w");
	FILE *fpmass = fopen("mass.txt","w");
	FILE *fptime = fopen("tm.txt","w");

/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	m  = 1.0;
	n  = 1.0;
	T = 2.0*pie/m;
	N = 100;



/////////////////////////////////////////RK things/////////////////////////////////////////


	double im_K_psi[2][ex_s][N],ex_K_psi[2][im_s][N];


/////////////////////////////// Box & res. setting ///////////////////////////////////////////////

	box_len = 2.0;
	dx = box_len/(double(N-1));

////////////////////////////// Time & dt settings ////////////////////////////////////////////////
	
	t_start = 0.0;
	t_end = 2.0*T;
	t_steps = 2000000;
	dt  = (t_end-t_start)/((double)t_steps);
	printf("dt %lf\n",dt);
	

////////////////////////////// Psi variables  /////////////////////////////////////////////////////

	double psi[N][2],lap_psi[2*N],lambda;

	int spd =N/4;
	int nN = N+ 2*spd;

	
	double k_grid[N+2*spd];
	
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

	initialise(psi,k_grid,dx,N,spd);

	int i,j;
	
	 	for(i=0;i<N;++i)
		{	
			
			fpGpsi[i][0] = psi[i][0] ;
			fpGpsi[i][1] = psi[i][1] ;

			if(i==(N-1))
		printf("iIni N-1 %lf %d\n",fpGpsi[i][0],i);
					
			
		

		}


///////////////////////  Evolution ///////////////////////////////////////////////////////////////

	int s_cntr,tcntr,printcntr,fail=0;
	printcntr = (int)(((double) t_steps)/100.0);	printf("%d\n",printcntr);
	double t,vel_val[2],c_psi[2],c_psi_amp,Vval,amp,avg_amp,amp_ini;
	double drc = 2.0*pie*n*2.0*pie*n;
	double fdt;



		
	

	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(1);t+=dt,++tcntr)
	{	
		
		avg_amp = 0.0;
		/*simple_lap(lap_psi,fpGpsi, dx,N);
		for(i=0;i<N;++i)
		{	dbi = (double)i;
			sol[0] = cos(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			sol[1] = -sin(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			

			//fdt = (fpGpsi[i-1][0]+fpGpsi[i+1][0]-2.0*fpGpsi[i][0])/(dx*dx);
		fprintf(fplap,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,lap_psi[2*i],lap_psi[2*i+1],K[i][1],-drc*fpGpsi[i][0],-drc*fpGpsi[i][1],0.0);
			
		}
		*/
		  fftw_execute(plan_pois_f);

	 	for(i=0;i<(N);++i)
		{	dbi = (double)(i);
			
			if((tcntr%printcntr)==0)  
			{
			  
			  amp  = sqrt(psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1]);				
			  avg_amp+=amp;
			  sol[0] = cos(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			  sol[1] = -sin(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			  fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp);
			  if(i==20)
			   fprintf(fptime,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp);
			}

			lambda = k_grid[i]*k_grid[i]*im_a[0][0]*dt/(2.0*m);


			fpGpsi_ft[i][0] = (fpGpsi_ft[i][0] + lambda*fpGpsi_ft[i][1])/(1.0+lambda*lambda);
			fpGpsi_ft[i][1] = (fpGpsi_ft[i][1] - lambda*fpGpsi_ft[i][0])/(1.0+lambda*lambda);

			fpGpsi_ft[i][0] = fpGpsi_ft[i][0]/((double)N);
			fpGpsi_ft[i][1] = fpGpsi_ft[i][1]/((double)N);

			/*K_ft[i][0] = -fpGpsi_ft[i][0]*(k_grid[i]*k_grid[i]);
			K_ft[i][1] = -fpGpsi_ft[i][1]*(k_grid[i]*k_grid[i]);

			if(i==0)
			{
			 K_ft[i][0] = 0.0;
			 K_ft[i][1] = 0.0;
			}
			*/
		}

		if(t==t_start)
		amp_ini = avg_amp;

		fftw_execute(plan_pois_b);
		//fftw_execute(plan_imp_b);
		simple_lap(lap_psi,fpGpsi, dx,N);
		

	

		for(s_cntr=1;s_cntr<imex_s;++s_cntr)
		{

			for(j=0;j<s_cntr;++j)
			{
			   for(i=0;i<N;++i)
			   {	
				
				
				if(j==0)
				{	
					if(s_cntr==1)
					{
					  c_psi[0] = psi[i][0];
		    			  c_psi[1] =  psi[i][1];
		

					}
					else
					{
					  c_psi[0] = fpGpsi[i][0];
		    			  c_psi[1] =  fpGpsi[i][1];
		

					}
	
		    			c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
		    			Vval = V(c_psi_amp);
		    			ex_vel(vel_val,c_psi,Vval);

		    			ex_K_psi[0][s_cntr-1][i] = vel_val[0];  ex_K_psi[1][s_cntr-1][i] = vel_val[1];


					

					im_K_psi[0][s_cntr-1][i] = -lap_psi[i*2+1]/(2.0*m);  im_K_psi[1][s_cntr-1][i] =  lap_psi[i*2]/(2.0*m);

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
				lambda = k_grid[i]*k_grid[i]*im_a[s_cntr][s_cntr]*dt/(2.0*m);

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
			 /// fftw_execute(plan_imp_b);
			simple_lap(lap_psi,fpGpsi, dx,N);

			
				
			


		}//RK stages concluded

		if((tcntr%printcntr)==0)  
		{
			  
			   printf("%lf\t%lf\n",t/t_end,avg_amp/((double)(N)));
		}
			
		avg_amp=0.0;

		for(i=0;i<N;++i)
		{	
	
			c_psi[0] = fpGpsi[i][0];
		    	c_psi[1] =  fpGpsi[i][1]; 
					
			c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
		    	Vval = V(c_psi_amp);
		    	ex_vel(vel_val,c_psi,Vval);

		    	ex_K_psi[0][imex_s-1][i] = vel_val[0];  ex_K_psi[1][imex_s-1][i] = vel_val[1];
			im_K_psi[0][imex_s-1][i] = -lap_psi[i*2+1]/(2.0*m);  im_K_psi[1][imex_s-1][i] =  lap_psi[i*2]/(2.0*m);
	
							
			   for(j=0;j<imex_s;++j)
			   {	psi[i][0]+=  dt*ex_b[j]*ex_K_psi[0][j][i]+ dt*im_b[j]*im_K_psi[0][j][i];
				psi[i][1]+=  dt*ex_b[j]*ex_K_psi[1][j][i]+ dt*im_b[j]*im_K_psi[1][j][i];
				
				

				
			    }


			fpGpsi[i][0] = psi[i][0];
			fpGpsi[i][1] = psi[i][1];


			


			amp  = sqrt(psi[i][0]*psi[i][0] + psi[i][1]*psi[i][1]);
			avg_amp+=amp;

			if((tcntr%printcntr)==0)  
			{ 
			  sol[0] = cos(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			  sol[1] = -sin(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			  fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,fpGpsi[i][0],fpGpsi[i][1],sol[0],sol[1],amp);
			}
			

		   }



		if((tcntr%printcntr)==0)  
		{ fprintf(fp,"\n\n\n");
		  fprintf(fp2,"\n\n\n");
		  fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)(N)));
		  printf("%lf\t%lf\n\n",t/t_end,avg_amp/((double)(N)));
		}
		


		

	}///// ENd f Time Evolution /////////////////////////
	avg_amp = 0.0;
	for(i=0;i<N;++i)
	{	amp  = sqrt(fpGpsi[i][0]*fpGpsi[i][0] + fpGpsi[i][1]*fpGpsi[i][1]);
		avg_amp+=amp;

	}

	fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)N));

	printf("Error in conserv. %lf\n", 100.0*fabs(avg_amp-amp_ini)/amp_ini);





	return(0);

}
