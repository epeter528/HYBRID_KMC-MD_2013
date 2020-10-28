#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include "kmc.h"


//int read_xvg_files(char *xvgfile,int *x,int sizeof_x,float *y,int sizeof_y)
int read_xvg_files(char *xvgfile,int *sizeof_x,float x[100000],float y[100000])
{
  
  int i,a,u;
  
  float step;
  
  float value;
  
  char char1 [200];
 
    
  
       FILE *fp = fopen(xvgfile,"r");
       
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
	
	*sizeof_x = u-1;	
	
	return 0; 

}

int read_xvg_files_int(char *xvgfile,int *sizeof_x,int x[500000],float y[500000])
{
  
  int i,a,u;
  
  int step;
  
  float value;
  
  char char1 [200];
  
       FILE *fp = fopen(xvgfile,"r");
       
       u = 1;
       	
        while(fgets(char1,sizeof(char1),fp) != 0) 
	{

	  if(memchr(char1,'#',strlen(char1)) !=NULL || memchr(char1,'@',strlen(char1)) !=NULL) {
	    
	//     printf("%s",char1);
	    
	  } else {	         
	    
	     sscanf(char1,"%d%f",&step,&value);
	 
	     x[u] = step;
	     
	     y[u] = value;
	    
       //     printf("%d\t%s\t%f\t%f\n",u,"u",step,value);	 
 
	     u ++;
	     
	  };	 
	  
	};
	
	*sizeof_x = u-1;	    
  
 
	
	return 0; 

}
