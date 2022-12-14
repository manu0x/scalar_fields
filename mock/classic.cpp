using namespace std;

#include <stdio.h>
#include <math.h>

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

	v[0] = 2.0*pie*(Vval*psi[1] - lap_psi[1]/(2.0*m*m))/T;
	v[1] = 2.0*pie*(-Vval*psi[0] + lap_psi[0]/(2.0*m*m))/T;



}


void initialise(double *psi,double dx,int N)
{

	int i,j; double di;
	for(i=0;i<N;++i)
	{
	   di = (double)i;	
	   *(psi+2*i) = sin(2.0*pie*di*n*dx);
	   *(psi+2*i+1) = 0.0;

	}


}


void lap(double *lpsi,double *psi, double dx,int N)
{
	int l1,l2,r1,r2,i;

	double vl1,vl2,vr1,vr2,vc;
	int Nm1 = N-1;
   for(i=0;i<(N-1);++i)
   {	l1 = (Nm1+ (i-1))%Nm1;
	l2 = (Nm1+ (i-2))%Nm1;

	r1 = (i+1)%Nm1;
	r2 = (i+2)%Nm1;


	vl1 = *(psi+2*l1);
	vl2 = *(psi+2*l2);

	vr1 = *(psi+2*r1);
	vr2 = *(psi+2*r2);

	vc = *(psi+2*i);
	*(lpsi+2*i) =  (-vr2 + 16.0*vr1 - 30.0*vc + 16.0*vl1 - vl2)/(12.0*dx*dx);


	vl1 = *(psi+2*l1+1);
	vl2 = *(psi+2*l2+1);

	vr1 = *(psi+2*r1+1);
	vr2 = *(psi+2*r2+1);

	vc = *(psi+2*i+1);
	*(lpsi+2*i+1) = (-vr2 + 16.0*vr1 - 30.0*vc + 16.0*vl1 - vl2)/(12.0*dx*dx);

    }		
	*(lpsi+2*i) = *(lpsi);
	*(lpsi+2*i+1) = *(lpsi+1);

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
	t_steps = 100000;
	dt  = (t_end-t_start)/((double)t_steps);
	printf("dt %lf\n",dt);
	
	

////////////////////////////// Psi variables  /////////////////////////////////////////////////////

	double psi[2*N],psik[5][2*N],lap_val[2*N];

	initialise(psi,dx,N);
///////////////////////  Evolution ///////////////////////////////////////////////////////////////

	int i,j,tcntr,printcntr,fail=0;
	printcntr = (int)(((double) t_steps)/10.0);	printf("%d\n",printcntr);
	double t,vel_val[2],c_psi[2],c_psi_amp,Vval,amp,avg_amp;

	for(i=0;i<N;++i)
	{
		*(psik[0]+2*i) = *(psi+2*i);
		*(psik[0]+2*i+1) = *(psi+2*i+1);
		
		

	}

		
	

	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(1);t+=dt,++tcntr)
	{	
		//if(tcntr%10==0) 
		//printf("time %lf %d\n",t/t_end,tcntr);

		for(j=1;j<=4;++j)
		{	
			
			lap(lap_val,&psik[0][0],dx,N);
			if(j==1)
			avg_amp = 0.0;
			for(i=0;i<N;++i)
			{	



			   if((tcntr%printcntr)==0 && (j==1) &&(!fail))  
			     {
				   
				   amp  = sqrt((*(psik[0]+2*i))*(*(psik[0]+2*i)) + (*(psik[0]+2*i+1))*(*(psik[0]+2*i+1)));
				   avg_amp+=amp;
				   fprintf(fp,"%lf\t%lf\t%lf\n",*(psi+2*i),*(psi+2*i+1),amp);
				   fprintf(fplap,"%lf\t%lf\t%lf\t%lf\n",i*dx,lap_val[2*i],lap_val[2*i+1],-(2.0*pie*n)*(2.0*n*pie)*sin(2.0*pie*n*i*dx));
				   
			     }
				

				
				c_psi[0] = *(psik[0]+2*i);
				c_psi[1] = *(psik[0]+2*i+1); 
				c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);
				Vval = V(c_psi_amp);
				vel(vel_val,c_psi,Vval,&lap_val[2*i]);
				
				*(psik[j]+2*i) = dt*vel_val[0];
				*(psik[j]+2*i+1) = dt*vel_val[1];
				if(j!=4)
				{ *(psik[0]+2*i) =  *(psi+2*i)+rk[j-1]*dt*vel_val[0];
				  *(psik[0]+2*i+1) = *(psi+2*i+1)+rk[j-1]*dt*vel_val[1];
				}

				else
				{ 
				  *(psi+2*i) = *(psi+2*i) + (1.0/6.0)*( *(psik[1]+2*i) + 2.0*(*(psik[2]+2*i)) + 2.0*(*(psik[3]+2*i)) + (*(psik[4]+2*i)) ); 
				  *(psi+2*i+1) = *(psi+2*i+1) + (1.0/6.0)*( *(psik[1]+2*i+1) + 2.0*(*(psik[2]+2*i+1)) + 2.0*(*(psik[3]+2*i+1)) + (*(psik[4]+2*i+1)) ); 

				  *(psik[0]+2*i) =  *(psi+2*i);
				  *(psik[0]+2*i+1) = *(psi+2*i+1);
				   amp  = sqrt((*(psik[0]+2*i))*(*(psik[0]+2*i)) + (*(psik[0]+2*i+1))*(*(psik[0]+2*i+1)));
				   if(isnan(amp))
				   {

					printf("%d failed at tcntr %d time %lf\n %lf %lf \n",i,tcntr,t/t_end,*(psik[0]+2*i),*(psik[0]+2*i+1));
					fail =1;
					break;
					
				   }	
				  
			            
				  
				}
			   
			}

			



		}

		if((tcntr%printcntr)==0)  
		{ fprintf(fp,"\n\n\n");
		  fprintf(fpmass,"%lf\t%lf\n",t/t_end,avg_amp/((double)(tcntr+1)));
		}
			            
				  

	}

	for(i=0;i<N;++i)
	{
		*(psik[0]+2*i) = *(psi+2*i);
		*(psik[0]+2*i+1) = *(psi+2*i+1);

		 amp  = sqrt((*(psik[0]+2*i))*(*(psik[0]+2*i)) + (*(psik[0]+2*i+1))*(*(psik[0]+2*i+1)));
		//fprintf(fp,"%lf\t%lf\t%lf\n",*(psi+2*i),*(psi+2*i+1),amp);

	}

	

	





	return(0);

}
