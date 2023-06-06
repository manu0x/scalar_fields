

using namespace std;
#include <stdio.h>



class imex_table
{
    protected:
    

    public:
    double **im_a,**ex_a;
    double *im_b,*ex_b;
    double *im_c,*ex_c;
    int s;

    imex_table(int stages)
    {
        int  i;
        s= stages;
       im_a = new double *[s];
         ex_a = new double *[s];
    
        for(i=0;i<s;++i)
        {
            im_a[i] = new double[s];
            ex_a[i] = new double[s];


        }
        im_b = new double[s];
        ex_b = new double[s];


    }

    void print_table()
    {   int i,j;
        printf("\nPrinting Implict A\n");
        for(i=0;i<s;++i)
        {
            for ( j = 0; j < s; j++)
            {
                printf("%lf ",im_a[i][j]);
            }
            printf("\n");

        }

        printf("\nPrinting Explict A\n");
        for(i=0;i<s;++i)
        {
            for ( j = 0; j < s; j++)
            {
                printf("%lf ",ex_a[i][j]);
            }
            printf("\n");

        }
        
        printf("\nPrinting Implict B\n");
        for(i=0;i<s;++i)
        {
            
            printf("%lf ",im_b[i]);
        
            printf("\n");

        }

        printf("\nPrinting Explict B\n");
        for(i=0;i<s;++i)
        {
            
            printf("%lf ",ex_b[i]);
        
            printf("\n");

        }


    }

    void read_from_file(char * param_file_name)
    {
        FILE *fp_param = fopen(param_file_name,"r");
        //printf("Opened file \n");
        int j,i=0;
        double num_val;
        int max_len = 80;
	    char read_c_line[max_len];
	    string read_str,cur_str, read_key, read_val;

        while(fgets(read_c_line,max_len,fp_param)!=NULL)
        {  read_str = read_c_line; //printf("re str  %s  %d\n",read_str.c_str(),read_str.length());
           if(read_str.length()>1)
           {
            
            for(j=0;j<s;++j)
            { 
                cur_str = read_str.substr(0,read_str.find(","));
                read_str = read_str.substr(read_str.find(",")+1,string::npos );

                remove_if(cur_str.begin(),cur_str.end(),::isspace);

                num_val = stod(cur_str);
                //printf("%d %d %lf\n ",i,j,num_val);

                if(i<s)
                {

                    im_a[i][j] = num_val;

                }
                else if (i<(2*s))
                {
                    ex_a[i-s][j] = num_val;
                }

                else if (i<(2*s+1))
                {
                    im_b[j] = num_val; 
                }
                else
                    ex_b[j] = num_val; 
                

               
            }
             ++i;
           }


        }

    }

};


/*
int main()
{


  imex_table ch(3);
 char nn[10] = "test";
 
  ch.read_from_file(nn);
  ch.print_table();


}
*/