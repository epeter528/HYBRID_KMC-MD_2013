#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kmc.h"

int barint(float *free_energy)

//int main (int argc, char *argv[])

{
   
  int i,j,k,l,m,n,u;
  
  int sizeof_x;
  
  float x[10000],y[10000];
  
  float fee_kt;
  
  FILE *fp,*fp2;
 
  char char1[300];
  
  float a,b,c,d,e;
  
  float step;
  
  float value;  
  // preprocess dhdl -files
  
  fp = fopen("dhdl1.xvg","r");

  fp2 = fopen("1.xvg","w");
  
  u = 1;
       	
  fprintf(fp2,"%s\n","@ s0 legend \"dH/d\\xl\\f{} \\xl\\f{} 0\"");
  fprintf(fp2,"%s\n","@ s1 legend \"\\xD\\f{}H \\xl\\f{} 0\"");
  fprintf(fp2,"%s\n","@ s2 legend \"\\xD\\f{}H \\xl\\f{} 0.5\"");
  fprintf(fp2,"%s\n","@ s3 legend \"\\xD\\f{}H \\xl\\f{} 1\"");
 
        while(fgets(char1,sizeof(char1),fp) != 0) 
	{

	  if(memchr(char1,'#',strlen(char1)) !=NULL || memchr(char1,'@',strlen(char1)) !=NULL) {
	    
	//     printf("%s",char1);
	    
	  } else {	         
	    
	    
	     fprintf(fp2,"%s",char1);

       //     printf("%d\t%s\t%f\t%f\n",u,"u",step,value);	 
 
     
	  };	 
	  
	}; 

  fclose(fp);
  fclose(fp2);
  
  fp = fopen("dhdl2.xvg","r");

  fp2 = fopen("2.xvg","w");
  
  u = 1;
       	
  fprintf(fp2,"%s\n","@ s0 legend \"dH/d\\xl\f{} \\xl\\f{} 0.5\"");
  fprintf(fp2,"%s\n","@ s1 legend \"\\xD\\f{}H \\xl\\f{} 0\"");
  fprintf(fp2,"%s\n","@ s2 legend \"\\xD\\f{}H \\xl\\f{} 0.5\"");
  fprintf(fp2,"%s\n","@ s3 legend \"\\xD\\f{}H \\xl\\f{} 1\"");
 
        while(fgets(char1,sizeof(char1),fp) != 0) 
	{

	  if(memchr(char1,'#',strlen(char1)) !=NULL || memchr(char1,'@',strlen(char1)) !=NULL) {
	    
	//     printf("%s",char1);
	    
	  } else {	         
	    
	    
	     fprintf(fp2,"%s",char1);

       //     printf("%d\t%s\t%f\t%f\n",u,"u",step,value);	 
 
     
	  };	 
	  
	}; 

  fclose(fp);
  fclose(fp2);	
  
  fp = fopen("dhdl3.xvg","r");

  fp2 = fopen("3.xvg","w");
  
  u = 1;
       	
  fprintf(fp2,"%s\n","@ s0 legend \"dH/d\\xl\f{} \\xl\\f{} 1\"");
  fprintf(fp2,"%s\n","@ s1 legend \"\\xD\\f{}H \\xl\\f{} 0\"");
  fprintf(fp2,"%s\n","@ s2 legend \"\\xD\\f{}H \\xl\\f{} 0.5\"");
  fprintf(fp2,"%s\n","@ s3 legend \"\\xD\\f{}H \\xl\\f{} 1\"");
 
        while(fgets(char1,sizeof(char1),fp) != 0) 
	{

	  if(memchr(char1,'#',strlen(char1)) !=NULL || memchr(char1,'@',strlen(char1)) !=NULL) {
	    
	//     printf("%s",char1);
	    
	  } else {	         
	    
	    
	     fprintf(fp2,"%s",char1);

       //     printf("%d\t%s\t%f\t%f\n",u,"u",step,value);	 
 
     
	  };	 
	  
	}; 

  fclose(fp);
  fclose(fp2);	
    
  system("~/GRO453/bin/g_bar -f 1.xvg 2.xvg 3.xvg -o t.xvg -temp 300.15");

       fp = fopen("t.xvg","r");
       
       u = 1;
       	
        while(fgets(char1,sizeof(char1),fp) != 0) 
	{

	  if(memchr(char1,'#',strlen(char1)) !=NULL || memchr(char1,'@',strlen(char1)) !=NULL) {
	    
	//     printf("%s",char1);
	    
	  } else {	         
	    
	     sscanf(char1,"%f%f",&step,&value);
	 
	     x[u] = step;
	     
	     y[u] = value;
	    
       //     printf("%d\t%s\t%f\t%f\n",u,"u",step,value);	 
 
	     u ++;
	     
	  };	 
	  
	};
	
	sizeof_x = u-1;
	
	
  fee_kt = 0;
  
  for(i=1;i<=sizeof_x;i++){
    
      fee_kt = fee_kt + y[i]; // sum up from lambda = 0, 0.5 , 1
    
  };
  
  *free_energy =  2.479*fee_kt; // recal into kJ/mol
  
//  printf("%f\t%s\n",free_energy,"free energy");
 
  return 0;
  
}
