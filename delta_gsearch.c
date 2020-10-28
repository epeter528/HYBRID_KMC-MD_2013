#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kmc.h"

int delta_g_search(int icounter,int scansteps,float energy_scan2[1000][7][50],int *execsteps,int *selectamino,
		   int *selectevent,float sim_temp,float OMEGA_B[7],float OMEGA_ZERO[7],float *delta_time,int number_of_event)

{
  int i,j,k,m,p,n;

  float PI = 3.141592654; 

  int maxbereich;
  
  int scanmin;
  
  int scanmax;
  
  float firstmax;
  
  float secondmin;
  
  float sumrate;
  
  int selectionnumber[1000];
  
  int selectionaminoacid[1000];
  
  int secondminscan;
  
  float firstmaximum[1000][7];
  
  int firstmax_scan[1000][7],sec_min_scan[1000][7];
  float min_step_en[1000][7];
  
  int exec_scanstep[1000][7];
  
  float firstmini[50];
  
  float rate_g[1000][7];
  float delta_g[1000][7];
  float secondminimum[1000][7];
   
  int  scanmini[50];
   
  int  endrate;
  
  float ran;
  
  float selectionrate[100]; 
  
  float friction_coeff = 1; 

  FILE *fp;

  
  for(i=1;i<=icounter;i++){
  
    *selectevent = 1;
    scanmin     = 1;
    
    for(j=1;j<=number_of_event;j++) {
    
 
       firstmax           =  energy_scan2[i][j][1];
       firstmaximum[i][j] =  energy_scan2[i][j][1];
       scanmax            = 1;

             for(k=1;k<=scansteps-1;k++){

               m = 0;
               maxbereich = 0; 

                 if(energy_scan2[i][j][k] > firstmax && energy_scan2[i][j][k+1] < energy_scan2[i][j][k]) {
		   
                 firstmax     = energy_scan2[i][j][k];
                 scanmax      = k;
                 m            = k; 
                 
		};
		
	     }; 

     firstmaximum[i][j]  = firstmax;
     firstmax_scan[i][j] = scanmax; 

      
    };

  };

  for(i=1;i<=icounter;i++){
    
    for(j=1;j<=number_of_event;j++){
      
      secondmin = energy_scan2[i][j][firstmax_scan[i][j]];
      secondminscan = 1;
      
          for(k=firstmax_scan[i][j];k<=scansteps-1;k++){
	
        	if(energy_scan2[i][j][k] < secondmin && energy_scan2[i][j][k+1] > energy_scan2[i][j][k]){
	  
	       secondmin = energy_scan2[i][j][k];
	       secondminscan = k;
	  
        	};
	
         };
    
         secondminimum[i][j] = secondmin;
         sec_min_scan[i][j]  = secondminscan;
      
         
         if(sec_min_scan[i][j] < firstmax_scan[i][j]){

         sec_min_scan[i][j] = firstmax_scan[i][j]; 

         };
	 
        };

  
        for(j=1;j<=number_of_event;j++){
  
	   min_step_en[i][j]   = energy_scan2[i][j][1];
          exec_scanstep[i][j] = sec_min_scan[i][j];


	  if(exec_scanstep[i][j] == 0 && firstmax_scan[i][j] == 0){
	    
	    firstmini[1] = energy_scan2[i][j][1];
	    scanmini[1]  = 1;
	    
	          p = 1;
		  m = 0;
		  
		  for(k=1;k<=scansteps;k++){
		    
		    if(k > 2 && energy_scan2[i][j][k-1] > energy_scan2[i][j][k] && energy_scan2[i][j][k] < energy_scan2[i][j][k+1] 
		      && energy_scan2[i][j][1] > energy_scan2[i][j][k]){
			
		      firstmini[p] = energy_scan2[i][j][k];
		      scanmini[p]  = k;
		    
		      p ++;
		    
		      m = k;
		      
		      
		      };
		    
		  };
		  
		  firstmax_scan[i][j] = scanmini[1];
		  firstmaximum[i][j]  = firstmini[1];
		  
		  for(k=m;k<=scansteps;k++){
		    
		    if(k > firstmax_scan[i][j] && energy_scan2[i][j][k] > firstmaximum[i][j]){
		      
		      firstmax_scan[i][j] = k;
		      firstmaximum[i][j]  = energy_scan2[i][j][k];

		      
		      n = k;
		      
		    };
		    
		    
		  };
		  
		  for(k=n;k<=scansteps;k++){
		    
		    exec_scanstep[i][j] = n;
		    secondminimum[i][j] = n;
		    
		      if(energy_scan2[i][j][k] < secondminimum[i][j]){
			
			secondminimum[i][j] = energy_scan2[i][j][k];
			exec_scanstep[i][j] = k;
			
			
		      };

		  };
	    
	  };
	  
        };       

	for(j=1;j<=number_of_event;j++){
	  
	  if(firstmaximum[i][j] == 0){
	    
	    firstmaximum[i][j] = energy_scan2[i][j][scansteps];
	    exec_scanstep[i][j] = scansteps;
	    
	  };
	  
	};
	
  };	

  for(i=1;i<=icounter;i++){
    
    for(j=1;j<=number_of_event;j++){
      
      rate_g[i][j] = 0;
      
    };
    
  };

  sumrate = 0;

fp = fopen("OUTPUT","a");
  
  for(i=1;i<=icounter;i++){
    
          fprintf(fp,"%d\t%s\n",i,"aminoacidnum");

    for(j=1;j<=number_of_event;j++){
      
      delta_g[i][j] = firstmaximum[i][j] - energy_scan2[i][j][1];
      
    //  if(delta_g[i][j] != 0) {
	
    //  rate_g[i][j] = (sqrt(pow(friction_coeff,2)/4 + pow(OMEGA_B[j],2)) - friction_coeff/2) * OMEGA_ZERO[j]/(2*PI*OMEGA_B[j])*exp(-delta_g[i][j]/(8.314*sim_temp));
          if(OMEGA_B[j] == 0) OMEGA_B[j] = 1;
      
	  rate_g[i][j] = (sqrt(pow(friction_coeff,2)/4 + pow(OMEGA_B[j],2)) - friction_coeff/2) * OMEGA_ZERO[j]/(2*PI*OMEGA_B[j])*exp(-delta_g[i][j]*6.022 /(8.314*sim_temp)); //energy diff in kT !	  

          if(isnanf(rate_g[i][j]) == 1) rate_g[i][j] = 1; 
	  
	  fprintf(fp,"%f\t%e\t%s\t%e\t%e\t%s\n" , delta_g[i][j] , rate_g[i][j] , "rate_g" , OMEGA_B[j] , OMEGA_ZERO[j] , "OMEGA B, OMEGA ZERO");
	  
   //   } else {
	
   //	rate_g[i][j] = 0;
	
   //   };
      
      sumrate = sumrate + rate_g[i][j];
      
    };
    
  };

fclose(fp);

k = 2;

selectionrate[1] = 0;

for(i=1;i<=icounter;i++) {
  
  for(j=1;j<=number_of_event;j++){
    
         if(isnanf(rate_g[i][j]) != 1){
    
              selectionrate[k] = selectionrate[k-1] + rate_g[i][j];
	      
	      selectionnumber[k] = j;
	      
	      selectionaminoacid[k] = i;
	      
	      k++;
	      
	 };
    
  };
  
};

endrate = k;

i = random_number(&ran);

for(i=1;i<=endrate;i++){
  
  if(selectionrate[i-1] < sumrate*ran && sumrate*ran < selectionrate[i]){
    
    *selectevent = selectionnumber[i];
    
    *selectamino = selectionaminoacid[i];
    
  };
  
};
 
     *execsteps = exec_scanstep[*selectamino][*selectevent];

i = random_number(&ran);

 //   sumrate = sumrate * 1E+15;

if(sumrate != 0) {

    *delta_time = - log(ran)/sumrate;
 
                 };
		 
if(sumrate == 0) *delta_time = 0;		 
    
//printf("%s\n","end of delta-g");
    
return 1;

}
