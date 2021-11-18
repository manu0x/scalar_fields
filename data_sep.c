#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define n 64

void main(int argc,char *args[])
{
	int i = atoi(args[1]);
	int tN =n*n*n;
	int j,k;
	printf("i is %d\n",i);

	FILE *fpin = fopen("psi.txt","r");

	char i_str[3];
	char f_out_str[7]="";
	sprintf(i_str,"%d",i);
	strcat(f_out_str,i_str); 
	strcat(f_out_str,".txt"); 
	
	//FILE *fp = fopen(".txt","w");
	printf("%s",f_out_str);
	
	FILE *fpout = fopen(f_out_str,"w");
	double a,dx0,dx1,dx2,psi_r_val,psi_i_val,dc;


	for(j=0;j<=((i+1)*tN);++j)
	{

		 fscanf(fpin,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&a,&dx0,&dx1,&dx2,&psi_r_val,&psi_i_val,&dc);

		if(j>=(i*tN))

		fprintf(fpout,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",a,dx0,dx1,dx2,psi_r_val,psi_i_val,dc);
	}
	


}
