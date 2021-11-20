#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void main(int argc,char *args[])
{
	int i = atoi(args[1]);
	int n = atoi(args[2]);
	int tN =n*n*n;
	int j,k;
	printf("i is %d %d\n",i,n);

	FILE *fpin = fopen("psi.txt","r");

	char i_str[3];
	char f_out_str[7]="";
	sprintf(i_str,"%d",i);
	strcat(f_out_str,i_str); 
	strcat(f_out_str,".txt"); 
	
	//FILE *fp = fopen(".txt","w");
	printf("%s",f_out_str);
	
	FILE *fpout = fopen(f_out_str,"w");
	double a,dx0,dx1,dx2,psi_r_val,psi_i_val,dc,psiamp2,a3;
	int ii=0;

	for(j=0;j<((i+1)*tN);++j)
	{

		 fscanf(fpin,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&a,&dx0,&dx1,&dx2,&psi_r_val,&psi_i_val,&psiamp2,&dc,&a3);

		if(j%(tN)==0)
		{	
			printf("ii %d\n",ii); ++ii;
		}

		if(j>=(i*tN))

		fprintf(fpout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",a,dx0,dx1,dx2,psi_r_val,psi_i_val,psiamp2,dc,a3);
	}
	


}
