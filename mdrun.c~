#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kmc.h"

int main(int argc, char *argv[])
{
  
   int seed;
  
   float x;
   
   seed = time(NULL);
   srand(seed);  
        
   srand( (unsigned)time( NULL ) );
   
   kmc(argc,argv);
  
return 0;
	 	 
}	 
int kmc(int argc, char *argv[])
{
         int mdtime,i,j,n,l,k,o,select_term,b,p;
	 
	 int range_1or2, buf[100],size;
	 
	 float Temp;

	 char char1[200];	 
	 
	 float T1,T2,T3,T4,T5,Temp_orig;
	 
	 float sim_temp;
	 
	 float  en_aver;
	 
	 float  noe_min;
	 
	 int sizeof_stddev;
	 
	 float stddev_dih[10000],rot_change;
	 
	 char *grofile,*topfile;
	 
	 char *trajfile,*final_file,*ener_file,*mdpfile;
	 
	 int aminoresnum;
	 
	 int scansteps;
	 
	 float energy_scan2[1000][7][50];
	 
	 int num_list1;
	 
	 float expconst;
	 
	 float delta_phi;
	 
	 int imp;

	 int icounter;
	 
	 int execsteps;
	 
	 int selectamino;
	 
	 int selectevent;
	 
	 int event_or_scan;
	 
	 int num_na,num_cl;

	 int replexbool;

	 int step;
	 
	 int kmcsteps;
	 
	 char chari[250];
	 
	 float sumrate;
	 
	 float OMEGA_B[7],OMEGA_ZERO[7];
	 
	 float omega_zero_val,num_eve_f;
	 
	 float delta_time;

         float scansteps2;
	 
         float restart_time; 

	 float friction_coeff,total_time;
	 
	 FILE *TIME;
	 
	 FILE *outkmc;
	 
	 FILE *MDTIMEFILE;
	 
	 Temp = 335.15;
	 
         int x[10000];
	 
	 int mdtime_max;
	 
	 int sum_of_mdtime;
      
         float y[10000],icount_f;
	 
	 char *xvgfile;
	 
	 float ran;
	 
	 int freezemax;

         int deppmd;
	 
	 int count_kmc;
	
         int number_of_events;
 
	 int range_checker_A, range_checker_B;

	 float eigenval2_x,eigenval2_y,eigenval2_z;
	 
	 float eig_x,eig_y,eig_z;
	 
	 float eig_x_array[100],eig_y_array[100],eig_z_array[100];
	
         float eig_x_array2[100],eig_y_array2[100],eig_z_array2[100];
 
	 int acc, acceptor[100],eventmax;

         float rota[100][10];
	 
	 int acceptor2,restart,restart_step;

         int cleanint;
	 
	 int evemax_int[100];

         int number_of_drystep;

	 FILE *fp = fopen("input","r");

	 //	 fscanf(fp,"%s",&chari);
	 //2
	 fscanf(fp,"%d",&mdtime_max);
	 
	 fscanf(fp,"%d",&freezemax);
         
	 
	 //3
//	 fscanf(fp,"%s",&chari);
         //4
	 fscanf(fp,"%d",&kmcsteps);
         //5
//	 fscanf(fp,"%s",&chari);
	 //5+1
	 fscanf(fp,"%d",&scansteps);
         //6
	 fscanf(fp,"%d",&replexbool);
         //7
//	 fscanf(fp,"%s",&chari);	 
         //8
	 fscanf(fp,"%d",&range_checker_A);
         //9
//	 fscanf(fp,"%s",&chari);
         //10
	 fscanf(fp,"%d",&range_checker_B);	 
         //11
//	 fscanf(fp,"%s",&chari);
         //12
	 fscanf(fp,"%d",&restart);	 
         //13
//	 fscanf(fp,"%s",&chari);
         //14
	 fscanf(fp,"%d",&cleanint);
         //15
//	 fscanf(fp,"%s",&chari);
         //16
	 fscanf(fp,"%d",&number_of_drystep);
         //17
//	 fscanf(fp,"%s",&chari);
         //18
	 fscanf(fp,"%f",&Temp_orig);
         //19
//	 fscanf(fp,"%s",&chari);
         //20
	 fscanf(fp,"%f",&noe_min);
         //21
//	 fscanf(fp,"%s",&chari);
         //22
	 fscanf(fp,"%d",&imp);
	 
	 fclose(fp);
	 
//	 printf("%d\n",restart);
	 
	 if(restart == 0) {
	 
	 TIME = fopen("TIME","w");
	 outkmc  = fopen("OUTPUT","w");
	 MDTIMEFILE = fopen("MDTIMEFILE","w");
	 
	 fprintf(MDTIMEFILE,"%s\n","start");
	 
	 fclose(MDTIMEFILE);
	 
//	 fscanf(fp,"%s",&chari);
	 //2
	 fprintf(outkmc,"%d\t%s\n",mdtime_max,"mdtime_max");
        //3
	 fprintf(outkmc,"%d\t%s\n",freezemax,"freezemax");
	 //	 fscanf(fp,"%s\n",&chari);
         //4
	 fprintf(outkmc,"%d\t%s\n",kmcsteps,"kmcsteps");
         //5
//	 fprintf(outkmc,"%s\n",&chari);
	 //5+1
	 fprintf(outkmc,"%d\t%s\n",scansteps,"scansteps");
         //6
	 fprintf(outkmc,"%d\t%s\n",replexbool,"replexbool");
         //7
//	 fprintf(outkmc,"%s\n",&chari);	 
         //8
	 fprintf(outkmc,"%d\t%s\n",range_checker_A,"range_checker_A");
         //9
//	 fprintf(outkmc,"%s\n",&chari);
         //10
	 fprintf(outkmc,"%d\t%s\n",range_checker_B,"range_checker_B");	 
         //11
//	 fprintf(outkmc,"%s\n",&chari);
         //12
	 fprintf(outkmc,"%d\t%s\n",restart,"restart");	 
         //13
//	 fprintf(outkmc,"%s\n",&chari);
         //14
	 fprintf(outkmc,"%d\t%s\n",cleanint,"cleanint");
         //15
//	 fprintf(outkmc,"%s\n",&chari);
         //16
	 fprintf(outkmc,"%d\t%s\n",number_of_drystep,"number_of_drystep");
         //17
//	 fprintf(outkmc,"%s\n",&chari);
         //18
	 fprintf(outkmc,"%f\t%s\n",Temp_orig,"Temporig");
         //19
//	 fprintf(outkmc,"%s\n",&chari);
         //20
	 fprintf(outkmc,"%f\t%s\n",noe_min,"noe_min");
         //21
//	 fprintf(outkmc,"%s\n",&chari);
         //22
	 fprintf(outkmc,"%d\t%s\n",imp,"imp");
	 
	 fclose(outkmc);
	 fclose(TIME);
	 
	 };
	
fp  = fopen ("mod.sh","w") ;

     fprintf(fp,"%s\n","#!/bin/bash");
     fprintf(fp,"%s\n","~/GRO453/bin/make_ndx -f minimized_water.gro -o index.ndx << EOF");
     fprintf(fp,"%s\n","q");
     fprintf(fp,"%s\n","<< EOF");
 
fclose(fp);

system("chmod 744 mod.sh");

system("./mod.sh");

if(restart == 0)  {	 
                  system("cp minimized_water.gro start.gro");
		   system("cp LOV2_bak.top start.top");	   
                  };

		   
sum_of_mdtime = 0;		  
		  
n = 1;	  
o = 1;

if(restart == 0)  total_time = 0;
		  
if(restart == 1)  {


          fp = fopen("TIME","r");

                       while(fgets(char1,sizeof(char1),fp) != NULL){

                         sscanf(char1,"%d%e",&l,&restart_time);
                         n = n + 1;

                             if(restart_time != 0) {

                                total_time = restart_time;

                                    };

                        };

          fclose(fp);


         fp = fopen("MDTIMEFILE","r");
 	  
                       while(fgets(char1,sizeof(char1),fp) != NULL){

                         o++;

                         if(o > 2) {
	  	    
			 sscanf(char1,"%d",&deppmd);
			 
                         };

                         if(deppmd > 0){

                            mdtime = deppmd;

                          };

                         if(deppmd > 0){

			 sum_of_mdtime = sum_of_mdtime + deppmd;

                         };
                
	                };
			
	  fclose(fp);		
  
	  
                  };		  
		  
//scansteps  = 20;
scansteps2 = (float)scansteps;


if(restart == 0) restart_step  = 1;
if(restart == 1) restart_step  = n;
if(restart == 1) {

   system("cp restart.top start.top");
   system("cp restart.gro start.gro");

   fp = fopen("OUTPUT","a");

   fprintf(fp,"%s\n","restarted ----- !");
   fprintf(fp,"%d\t%s\n",restart_step,"step");
   fprintf(fp,"%d\t%s\n",sum_of_mdtime,"sumofmdtime");
   fprintf(fp,"%s\t%d\n","mdtime",mdtime);
   fprintf(fp,"%s\t%e\n","total_time",total_time);

   fclose(fp);

   system("~/GRO453/bin/grompp -f minim.mdp -c start.gro -p start.top -o t.tpr");
   system("~/GRO453/bin/mdrun  -s t.tpr -c start.gro -p start.top ");

   fp = fopen("OUTPUT","a");

   fprintf(fp,"%s\n","minimized");

   fclose(fp);

};


for(l = restart_step; l <= kmcsteps; l++){

    if(l == 1 && restart == 0) {
      
      mdtime = 10000;
      
    }; 
  

    if(l > 1) {

      mdtime = mdtime*1.5;

    };

    
    if(mdtime > mdtime_max){
      
       mdtime = mdtime_max;
      
    };
    
    i = freeze(freezemax);
    
 

    
    if(replexbool == 1){    
      
      
      i = replex(argc,argv,T1,T2,T3,T4,T5,&sim_temp,mdtime,Temp_orig,imp);

      
    };
      
    
    if(replexbool == 0){
      
      
      sim_temp = Temp_orig;
      
    };
    
    grofile = "start.gro";
   

    i = mdprod(argc,argv,mdtime,sim_temp,imp,grofile);


    MDTIMEFILE = fopen("MDTIMEFILE","a");

    fprintf(MDTIMEFILE,"%d\n",mdtime);

    fclose(MDTIMEFILE);

    
    outkmc = fopen("OUTPUT","a");

    fprintf(outkmc,"%s\n","End of MD");
  
    fclose(outkmc);
  
    if( l == 1 ) {    
  
    i = solvent_config_separator(&num_na,&num_cl,imp);         
    
    };

    i = edit(argc,argv); 

    i = noe_sep(argc,argv,range_checker_A,range_checker_B,eventmax);
    
    i = dih_trans(argc,argv,stddev_dih,&sizeof_stddev);
    
    xvgfile = "searchlist.ndx";  
  
    i = read_xvg_files_int(xvgfile,&icounter,x,y);
    
    outkmc = fopen("OUTPUT","a");   

    
    fprintf(outkmc,"%s\t%d\n","Begin of Trans-scan",icounter); 
    
    fclose(outkmc); 
    
    
    i = kramer_zero(argc,argv,&omega_zero_val,sim_temp,aminoresnum);

    for(i=1;i<=7;i++){

        if(omega_zero_val == 0) omega_zero_val = 1;

        OMEGA_ZERO[i] = omega_zero_val;
        

    };
    
    event_or_scan = 0; //kramer switched on !    

    i = random_number(&ran);
    
    if(icounter <= 0){
    
    icount_f = (float) icounter;

    count_kmc = roundf(icount_f*ran*2);
    
    if(count_kmc == 0) {
      
      random_number(&ran);

      count_kmc = roundf(ran*freezemax+1);
      
    };
    
    };
    
    if(icounter > 0){
      
       count_kmc = 1;
      
    };
      
    for(o = 1; o <= count_kmc; o++){ 
      
      
      if(o >= 2) system("cp minimized.gro minimized3.gro");

        i = random_number(&ran);  

 //       number_of_events = (int)roundf(ran*6)+1;

        number_of_events = 7; // previous ran number of events makes no sense

       //   number_of_events = 7;
       
        if(icounter > 4) icounter = 4; 

        for(k = 1; k <= icounter; k++)
      
             {     
      
               aminoresnum = x[k];
      
               if(stddev_dih[k] == 0 || isnanf(stddev_dih[k]) == 1) stddev_dih[k] = 2;
      

	       for(b = 1; b <= number_of_events; b ++) {	  
		  

                   i = random_number(&ran);

                 /*

                       rot_change  = stddev_dih[k] *ran ;

                       delta_phi   = stddev_dih[k] *ran ;


                  if(delta_phi*scansteps2 >= 2.5) {

                        delta_phi  = 2.5 * ran ;

                   };

                 if(rot_change*scansteps2 >= 2.5) {

                        rot_change = 2.5 * ran ;

                   }; */ // this is now switched of because the adaptive range of dih was maybe nonsense.

                  delta_phi = 0.0025;

                 rot_change = 0.0025;

                  rota[k][b] = rot_change ;
		 

		      num_eve_f = (float) number_of_events;

		//      i = random_number(&ran);
		 
		//      evemax_int[b] = (int)roundf(ran*6) ;

                //  This is changed here because random selection seems to make no sense !

                evemax_int[b] = b;

      // determination of middle point in distribution ... max of event, around of which we find events
		      
		// if !      
		/*      if(evemax_int[b] == 1)
                                                 {
			
			                          i = dihedralbreak(argc,argv,scansteps,aminoresnum,rot_change,energy_scan2,k,sim_temp,OMEGA_B,event_or_scan,b);      
      
                                                 outkmc = fopen("OUTPUT","a");      
	 
                                                    fprintf(outkmc,"%s\t%f\n","dih-scan",rot_change);
       
                                                 fclose(outkmc); 		     
		     
		                                  };  */   
	 
                 // if !
		      if(evemax_int[b] == 1)		    
		                                  {
		    
                                                 expconst = y[k]/3;		    

						   i = break_trans(argc,argv,scansteps,aminoresnum,energy_scan2,k,sim_temp,expconst,OMEGA_B,event_or_scan,
								   eigenval2_x,eigenval2_y,eigenval2_z,&eig_x,&eig_y, &eig_z,b);
		  
						   
                                                 outkmc = fopen("OUTPUT","a"); 
   
                                                    fprintf(outkmc,"%s\n","break-scan");
       
                                                 fclose(outkmc); 						   
		  
						  };
		  
	                                       	  eig_x_array[k] = eig_x;
		                                         eig_y_array[k] = eig_y;
	                                           	  eig_z_array[k] = eig_z;		  
              
       
		// if! 
		      if(evemax_int[b] == 2)							  
					        {
		 
                                                i = form_trans(argc,argv,scansteps,aminoresnum,energy_scan2,k,sim_temp,expconst,OMEGA_B,event_or_scan,
							       eigenval2_x,eigenval2_y,eigenval2_z,&eig_x,&eig_y,&eig_z,b);
		 
	
						      outkmc = fopen("OUTPUT","a");        
       
                                                            fprintf(outkmc,"%s\n","form-scan");  
       
                                                    fclose(outkmc);
						
						};
						
	                                         	eig_x_array2[k] = eig_x;
                                                      eig_y_array2[k] = eig_y;
                                                      eig_z_array2[k] = eig_z;
						      
               // if !
		      if(evemax_int[b] == 3)						      
						      
		                               {
       
                                               i = phi_trans(argc,argv,scansteps,aminoresnum,energy_scan2,k,delta_phi,sim_temp,OMEGA_B,event_or_scan,b);
       
				                outkmc = fopen("OUTPUT","a");        
        
                                                        fprintf(outkmc,"%s\n","phi-scan");  
       
                                               fclose(outkmc); 	       
					       
					       
					        };
			
		// if !				
		      if(evemax_int[b] == 4)						
	                                     {
	       
	       
                                              i = psi_trans(argc,argv,scansteps,aminoresnum,energy_scan2,k,delta_phi,sim_temp,OMEGA_B,event_or_scan,b);
       
                                              outkmc = fopen("OUTPUT","a");
		
                                                        fprintf(outkmc,"%s\n","psi-scan"); 
       
                                              fclose(outkmc);  					      
					      
					     };
					      
		// if !	
		      if(evemax_int[b] == 5)					     
       
	                                     {
	       
                                              i = phi_trans_neg(argc,argv,scansteps,aminoresnum,energy_scan2,k,delta_phi,sim_temp,OMEGA_B,event_or_scan,b); 
       
                                             outkmc = fopen("OUTPUT","a");        
        
                                                      fprintf(outkmc,"%s\n","phi-scan-neg");   
       
                                             fclose(outkmc); 					      
					      
					      
					     };        
       
		// if!	
		      if(evemax_int[b] == 6)
			
	                                    {
	       
                                            i = psi_trans_neg(argc,argv,scansteps,aminoresnum,energy_scan2,k,delta_phi,sim_temp,OMEGA_B,event_or_scan,b);       

					      
                                              outkmc = fopen("OUTPUT","a");       

                                                    fprintf(outkmc,"%s\n","psi-scan-neg");
       
                                              fclose(outkmc);  					    
					    
					     };
					     
		      if(evemax_int[b] == 7 && k == 1) // need to call rep-stretch only once ! we can copy the energies times icounter !
			
	                                    {
	       
                                            i = repeated_stretch(argc,argv,scansteps,range_checker_A,range_checker_B,energy_scan2,k,sim_temp,expconst,OMEGA_B,event_or_scan,
								   eigenval2_x,eigenval2_y,eigenval2_z,&eig_x,&eig_y, &eig_z,b);       
				      
                                              outkmc = fopen("OUTPUT","a");       

                                                    fprintf(outkmc,"%s\n","repeated_stretch");
       
                                              fclose(outkmc);  					    
					    
					      for(p = 2; p <= icounter; p ++){
						
						       for(o = 1; o <= scansteps; o++){
						
						          energy_scan2[p][7][o] = energy_scan2[1][7][o]; 
						
						         }; 
							 
					           };
					      
					      
					     };					     
					    
      
	       }; // end of loop over number_of_events loop_number = b;
	       
       }; // End of Transition-state sampling
  
      
   
       outkmc = fopen("OUTPUT","a");   

      fprintf(outkmc,"%s\n","End of trans-state sampling");

      fclose(outkmc);  
   
   
      i = delta_g_search(icounter,scansteps,energy_scan2,&execsteps,&selectamino,
		   &selectevent,sim_temp,OMEGA_B,OMEGA_ZERO,&delta_time,number_of_events);

     if(cleanint == 1 && l % 2 == 0){ 
 
        system("rm -f ~/.slurm/* ");
 
     };

     total_time = total_time + delta_time;
   
     printf("%e\n",total_time);
     
     TIME = fopen("TIME","a");
     
     fprintf(TIME,"%d\t%e\n",l,total_time/1E+15);

     fclose(TIME);
     
     outkmc = fopen("OUTPUT","a");      
     
     fprintf(outkmc,"%s\n","---");     
     fprintf(outkmc,"%s\n","DELTA-TIME");
     fprintf(outkmc,"%e\n",delta_time/1E+15);
     fprintf(outkmc,"%s\n","TOTAL-TIME");
     fprintf(outkmc,"%e\n",total_time/1E+15);
     fprintf(outkmc,"%s\n","select-event");
     fprintf(outkmc,"%d\n",selectevent);
     fprintf(outkmc,"%s\n","execsteps");
     fprintf(outkmc,"%d\n",execsteps);
     fprintf(outkmc,"%s\n","selectamino");
     fprintf(outkmc,"%d\n",x[selectamino]);
     fprintf(outkmc,"%s\t%d\n","step",l);
     fprintf(outkmc,"%s\t%d\n","kmc-int",o);
     fprintf(outkmc,"%s\t%d\n","kmc-cyc",count_kmc);     
     fprintf(outkmc,"%s\n","---");

     fclose(outkmc);   
     
     rot_change = rota[selectamino][selectevent];

     delta_phi  = rota[selectamino][selectevent]; 

     event_or_scan = 1;
     
//     if(evemax_int[selectevent] == 1){
      
//          i = dihedralbreak(argc,argv,execsteps,x[selectamino],rot_change,energy_scan2,1,sim_temp,OMEGA_B,event_or_scan,b);       
       
//    };
	   
     if(evemax_int[selectevent] == 1){
      
           eigenval2_x = eig_x_array[selectamino];
           eigenval2_y = eig_y_array[selectamino];
           eigenval2_z = eig_z_array[selectamino];
      
          i = break_trans(argc,argv,execsteps,x[selectamino],energy_scan2,1,sim_temp,expconst,OMEGA_B,event_or_scan,eigenval2_x,eigenval2_y,eigenval2_z,&eig_x,&eig_y, &eig_z,b);
      
    };
    
     if(evemax_int[selectevent] == 2){
      
           eigenval2_x = eig_x_array2[selectamino];
           eigenval2_y = eig_y_array2[selectamino];
           eigenval2_z = eig_z_array2[selectamino];

      
          i = form_trans(argc,argv,scansteps,aminoresnum,energy_scan2,k,sim_temp,expconst,OMEGA_B,event_or_scan,eigenval2_x,eigenval2_y,eigenval2_z,&eig_x,&eig_y,&eig_z,b);
      
    };
    
     if(evemax_int[selectevent] == 3){
      
          i = phi_trans(argc,argv,execsteps,x[selectamino],energy_scan2,1,delta_phi,sim_temp,OMEGA_B,event_or_scan,b); 
      
    };
    
     if(evemax_int[selectevent] == 4){
      
          i = psi_trans(argc,argv,execsteps,x[selectamino],energy_scan2,1,delta_phi,sim_temp,OMEGA_B,event_or_scan,b);

      
    };
    
     if(evemax_int[selectevent] == 5){
      
          i = phi_trans_neg(argc,argv,execsteps,x[selectamino],energy_scan2,1,delta_phi,sim_temp,OMEGA_B,event_or_scan,b);
      
    };
    
     if(evemax_int[selectevent] == 6){
      
          i = psi_trans_neg(argc,argv,execsteps,x[selectamino],energy_scan2,1,delta_phi,sim_temp,OMEGA_B,event_or_scan,b);
      
    };  
    
      if(evemax_int[selectevent] == 7){ // need to call rep-stretch only once ! we can copy the energies times icounter !
	       
          i = repeated_stretch(argc,argv,scansteps,range_checker_A,range_checker_B,energy_scan2,k,sim_temp,expconst,OMEGA_B,event_or_scan,
	                       eigenval2_x,eigenval2_y,eigenval2_z,&eig_x,&eig_y, &eig_z,b);     
	};
					    // end of KMC_scan

  }; //end of loop over all KMC-scans of one phase
    

     i = traj(argc,argv,sum_of_mdtime,l,imp);
     
   if( l < 0 ){ 
 
     i = solvent(argc,argv,imp);
        
   };

   if( l > number_of_drystep){

     i = solvent2(argc,argv,imp);

   };

    sum_of_mdtime = sum_of_mdtime + mdtime;
    
}; // End of KMC-MAINLOOP
	 	 
fclose(TIME);
fclose(outkmc);

	 
}

 int freeze(int freezemax){
   
   float x;
   
   int i,l;
   
   int n,nfreeze;
   
   int a;  
  
  char char1[100];
      
  char res_type[5],atom_type[2];
  
  int atom_num;
  
  float atom_x,atom_y,atom_z,freezemax_float;

  int res_number;  
  
  int natoms,atom_tot;
  
  FILE *fp;
  FILE *fp2;
  
  float box_x,box_y,box_z,natoms_fl;
  
  char *chari;
  
  int count_na,count_cl;

  int freezeatom[1000];
  
  int atomfreeze[1000];
  
  int number_of_freezeres;
 
      fp = fopen("start.gro","r");
      
      a = 0;
      
      natoms = 0;
      
      count_na = 0;
      
      count_cl = 0; 
      
      n = 1;
      
      while(fgets(char1,sizeof(char1),fp) != NULL){
	
	            a ++;
	
       //	            printf("%s\n",char1);
		    
		if(a ==2) {
		  
		sscanf(char1,"%d",&atom_tot);
		
                };
	
		if(a > 2 && a <= atom_tot+2) {    
	            sscanf(char1,"%5d%5s%2s%5d%7f%7f%7f",&res_number,res_type,atom_type,&atom_num,&atom_x,&atom_y,&atom_z);

	//	    printf("%s\n",atom_type);

                    if(strcmp(res_type,"TRP") ==0 || strcmp(res_type,"VAL") == 0 
                        ||     strcmp(res_type,"ILE") == 0 ||
                               strcmp(res_type,"LEU") == 0 ||
                               strcmp(res_type,"PHE") == 0 ||
                               strcmp(res_type,"CYS") == 0 ||
                               strcmp(res_type,"TYR") == 0) {

                      if(strcmp(atom_type,"C")  != 0 &&
                         strcmp(atom_type,"N")  != 0 &&
                         strcmp(atom_type,"CA") != 0 &&
                         strcmp(atom_type,"O")  != 0) {

                          atomfreeze[n] = atom_num ;

                          n++;

                        };
	   
                   }; 
	
                   if(strcmp(atom_type,"Na")==0){
	    
	            count_na++;
	     
	            };	            

                   if(strcmp(atom_type,"Cl")==0){
	    
	            count_cl++;
	    
	           };	
		   
		};	
		
		if(a == atom_tot+3) sscanf(char1,"%f%f%f",&box_x,&box_y,&box_z);
	
      };
      
      fclose(fp);   
      
      number_of_freezeres = n-1;
      
      fp = fopen("posre_spec.itp","w");
      
      fprintf(fp,"%s\n","[ position_restraints ]");
   
     i = random_number(&x);
   
     freezemax_float = (float)(freezemax);
     
     printf("%f\n",x);
     
     nfreeze = roundf(x*freezemax_float);
     
     printf("%d\t%s\n",nfreeze,"nfreeze");
     
     natoms_fl = (float) natoms;
     
     for(i=1;i<=nfreeze;i++){
       
         n = random_number(&x);
	 
	 freezeatom[i] = atomfreeze[(int)roundf(x*number_of_freezeres)];
       
    };
    
     for(i=1;i<=nfreeze-1;i++){
       
         for(n = 2;n <= nfreeze; n++){
	  
	   if(freezeatom[i] == freezeatom[n]){
	     
	      l = random_number(&x);  
	    
	      freezeatom[i] = freezeatom[i] + i;
	     
	  };
	   
	};
       
    };    
     
     for(i=1;i<=nfreeze;i++){
       
         fprintf(fp,"%d\t%s\n",freezeatom[i]," 1   50  50  50 ");
	 
	 printf("%d\n",freezeatom[i]);
       
    };
    
    fclose(fp);

    system("cp posre_spec.itp posre.itp");
    
}

 int edit(int argc,char *argv[])
 {
 
   int i,a;
   
/*   argc = 8 ;
   
   argv[0] = "editconf";
   
   argv[1] = "-f";
   
   argv[2] = "minimized3.gro";
   
   argv[3] = "-c";
   
   argv[4] = "-o";
   
   argv[5] = "minimized3.gro";
   
   argv[6] = "-resnr";
   
   argv[7] = "1"; */
   
  // i = gmx_editconf2(argc,argv,a); 
  // system("editconf -f minimized3.gro -c -o minimized3.gro -resnr 1");

   system("rm -f ./#*#");
   
}

 int traj(int argc,char *argv[],int mdtime,int step,int imp)
 
 {
   
   int i,a;
   
   int atom_tot;
   
   char char1[500];
   
   char res_type[5],atom_type[5];
  
   int atom_num;
  
   float atom_x,atom_y,atom_z;
   
   float mdtime_float,step_float;

   int res_number; 
   
   FILE *fp,*fp2;
   
   mdtime_float = (float)mdtime;
   step_float   = (float)step;
  
   mdtime_float = mdtime_float/1000;   
   
  if(imp == 0) {
    
    
      fp = fopen("minimized3.gro","r");
      fp2 = fopen("prot.ndx","w");
      
      fprintf(fp2,"%s\n","[ protein ]");
      
      a = 0;
      
      
      while(fgets(char1,sizeof(char1),fp) != NULL){
	
	            a ++;
	
       //	            printf("%s\n",char1);
		    
		if(a ==2) {
		  
		sscanf(char1,"%d",&atom_tot);
		
                };
	
		if(a > 2 && a <= atom_tot+2) {    
	            sscanf(char1,"%d%s%s%d%f%f%f",&res_number,res_type,atom_type,&atom_num,&atom_x,&atom_y,&atom_z);

	//	    printf("%s\n",atom_type);
	
                    if(strcmp(res_type,"SOL")!=0){
	    
	                 fprintf(fp2,"%d\n",atom_num);
	    
	            };

		  
		};	
	
      };
      
      fclose(fp); 
      fclose(fp2);
      
   //   printf("%d\n",count_na);
   //   printf("%d\n",count_cl);
      
   //   printf("%d\t%s\n",natoms,"natoms");
    
/*  argc = 7;
  argv[0] ="trjconv";
  argv[1] ="-f";
  argv[2] ="md_traj.trr";
  argv[3] ="-o";
  argv[4] ="md_traj.trr";
  argv[5] ="-n";
  argv[6] ="prot.ndx";*/
   
//  i = gmx_trjconv(argc,argv); 

// system("trjconv -s run.tpr -f md_traj.trr -o md_traj.trr -n prot.ndx");

  };
  
  if(step == 1) system("mv md_traj.trr md_trajx_internal.trr");
  
  if(step >= 2) {
  
 /*  argc = 7;
   
   argv[0] = "trjcat";
   argv[1] = "-f";
   argv[2] = "md_trajx_internal.trr";
   argv[3] = "md_traj.trr";
   argv[4] = "-settime";
   argv[5] = "-o";
   argv[6] = "md_trajx_internal.trr";
   
   i = gmx_trjcat2(argc,argv,step,mdtime);*/
 
   fp = fopen("trjcat.sh","w");
   
   fprintf(fp,"%s\n","#!/bin/bash");
   fprintf(fp,"%s\n","~/GRO453/bin/trjcat -f md_trajx_internal.trr md_traj.trr -settime -o md_trajx_internal.trr << EOF");
   fprintf(fp,"%s\n","0");
   fprintf(fp,"%f\n", mdtime_float);
   fprintf(fp,"%s\n","EOF");
   
   fclose(fp);
   
   system("chmod 744 ./trjcat.sh");
   
   system("./trjcat.sh");
   
   system("rm -f ./#*#");
   
  };
   
}

 int solvent2(int argc,char *argv[],int imp)
 
 {
   FILE *fp;
   FILE *fp2;

    fp2 = fopen("minim.mdp","w");
	
	fprintf(fp2,"%s\n","title                    = implicit minim.");
	fprintf(fp2,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp2,"%s\n","include 		 = -I../top");
        fprintf(fp2,"%s\n","define               = -DPOSRE");
        fprintf(fp2,"%s\n","integrator		 = l-bfgs");
	fprintf(fp2,"%s\n","emtol                = 8000");
	fprintf(fp2,"%s\n","emstep               = 0.005");
	fprintf(fp2,"%s\n","nsteps               = 200");
	fprintf(fp2,"%s\n","dt			 = 0.001");
	fprintf(fp2,"%s\n","comm_mode            =  linear ");
	fprintf(fp2,"%s\n","nstxout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstvout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp2,"%s\n","nstenergy		 = 1500 ");
	fprintf(fp2,"%s\n","nstxtcout		 = 1500 ");	
	fprintf(fp2,"%s\n","xtc_grps		 = protein");
	fprintf(fp2,"%s\n","energygrps		 = protein ");	
	fprintf(fp2,"%s\n","nstlist 		 = 10");	
	fprintf(fp2,"%s\n","ns_type 		 = grid");		
	fprintf(fp2,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp2,"%s\n","coulombtype		 = pme");		
	fprintf(fp2,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp2,"%s\n","vdwtype           = shift ");	
	fprintf(fp2,"%s\n","rvdw			 = 1.0");	
	fprintf(fp2,"%s\n","pbc                  = xyz");	
	fprintf(fp2,"%s\n","tcoupl  		 = no");	
	fprintf(fp2,"%s\n","tc-grps 		 = system");	
	fprintf(fp2,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp2,"%s\n","ref_t		 = 300");	
	fprintf(fp2,"%s\n","Pcoupl  		 = no");	
	fprintf(fp2,"%s\n","gen_vel 		 = no");	
	fprintf(fp2,"%s\n","gen_temp		 = 300");	
        fprintf(fp2,"%s\n","constraints         = none");
	fclose(fp2);	
  	 
    system("~/GRO453/bin/grompp -f minim.mdp -c minimized.gro -p start.top -o run.tpr");

    system("~/GRO453/bin/mdrun -s run.tpr -c start.gro ");

    system("cp ./start.gro ./restart.gro");
    system("cp ./start.top ./restart.top");
   
   return 1;

};

 int solvent(int argc,char *argv[],int imp)
 
 {
   FILE *fp;
   FILE *fp2;
   
   int a,i;
   
   char char1[500];
   
   float box_x,box_y,box_z;
   
   float max_x,max_y,max_z;
   
   float min_x,min_y,min_z;
   
   box_x = box_y = box_z = 10;
 
   char *grofile;
     
   char *topfile;
	
   char  *mdpfile; 
	
   char  *ndxfile;
   
   char res_type[5],atom_type[2];
  
   int atom_num;
  
   float atom_x,atom_y,atom_z;

   int res_number;  
  
   int natoms,atom_tot;
   
   char str[100];
 
  char char_s[500];  
   
/*   fp = fopen("minimized.gro","r");
      
      a = 0;
      
      natoms = 0;
      
      max_x = max_y = max_z = 0;
      
      while(fgets(char1,sizeof(char1),fp) != NULL){
	
	            a ++;
	
       //	            printf("%s\n",char1);
		    
		if(a ==2) {
		  
		sscanf(char1,"%d",&natoms);
		
                };  
		
		if(a >=2 && a <= natoms+2){
		  
		  sscanf(char1,"%d%s%s%d%f%f%f",&res_number,res_type,atom_type,&atom_num,&atom_x,&atom_y,&atom_z);
		  
		  if(atom_x > max_x) max_x = atom_x;
		  if(atom_y > max_y) max_y = atom_y;
		  if(atom_z > max_z) max_z = atom_z;		  
	//	  printf("%d\t%d\t%d\t%f\t%f\t%f\n",natoms,a,res_number,atom_x,atom_y,atom_z);
		  
		};
   
      };
      
   fclose(fp);   

   fp = fopen("minimized.gro","r");
      
      a = 0;
      
      natoms = 0;
      
      min_x = max_x;
      min_y = max_y;     
      min_z = max_z;     
      
      while(fgets(char1,sizeof(char1),fp) != NULL){
	
	            a ++;
	
       //	            printf("%s\n",char1);
		    
		if(a ==2) {
		  
		sscanf(char1,"%d",&natoms);
		
                };  
		
		if(a >=2 && a <= natoms+2){
		  
		  sscanf(char1,"%d%s%s%d%f%f%f",&res_number,res_type,atom_type,&atom_num,&atom_x,&atom_y,&atom_z);
		  
		  if(atom_x < min_x) min_x = atom_x;
		  if(atom_y < min_y) min_y = atom_y;
		  if(atom_z < min_z) min_z = atom_z;		  
	//	  printf("%d\t%d\t%d\t%f\t%f\t%f\n",natoms,a,res_number,atom_x,atom_y,atom_z);
		  
		};
   
      };
      
   fclose(fp);  */
   
//   sprintf(str,"%f\t%f\t%f\t",(max_x - min_x)+3.5,(max_y - min_y)+3.5,(max_z - min_z)+3.5); 

   system("~/GRO453/bin/editconf -f minimized.gro -o minimized.gro -c ");

   system("cp LOV2.top LOV2gen.top");
   
   if(imp == 0){
     

     system("~/GRO453/bin/genbox -cp minimized.gro -cs spc216.gro -p LOV2gen.top -o start.gro");
     
     
   };
   
   printf("%s\n","here after genbox");
  

   system("cp LOV2gen.top start.top");
   
 //  system("rm -f ./#*#");

    fp2 = fopen("minim.mdp","w");
	
	fprintf(fp2,"%s\n","title                    = implicit minim.");
	fprintf(fp2,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp2,"%s\n","include 		 = -I../top");
        fprintf(fp2,"%s\n","define               = -DPOSRE");
        fprintf(fp2,"%s\n","integrator		 = l-bfgs");
	fprintf(fp2,"%s\n","emtol                = 8000");
	fprintf(fp2,"%s\n","emstep               = 0.005");
	fprintf(fp2,"%s\n","nsteps               = 200");
	fprintf(fp2,"%s\n","dt			 = 0.001");
	fprintf(fp2,"%s\n","comm_mode            =  linear ");
	fprintf(fp2,"%s\n","nstxout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstvout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp2,"%s\n","nstenergy		 = 1500 ");
	fprintf(fp2,"%s\n","nstxtcout		 = 1500 ");	
	fprintf(fp2,"%s\n","xtc_grps		 = protein");
	fprintf(fp2,"%s\n","energygrps		 = protein ");	
	fprintf(fp2,"%s\n","nstlist 		 = 10");	
	fprintf(fp2,"%s\n","ns_type 		 = grid");		
	fprintf(fp2,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp2,"%s\n","coulombtype		 = pme");		
	fprintf(fp2,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp2,"%s\n","vdwtype           = shift ");	
	fprintf(fp2,"%s\n","rvdw			 = 1.0");	
	fprintf(fp2,"%s\n","pbc                  = xyz");	
	fprintf(fp2,"%s\n","tcoupl  		 = no");	
	fprintf(fp2,"%s\n","tc-grps 		 = system");	
	fprintf(fp2,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp2,"%s\n","ref_t		 = 300");	
	fprintf(fp2,"%s\n","Pcoupl  		 = no");	
	fprintf(fp2,"%s\n","gen_vel 		 = no");	
	fprintf(fp2,"%s\n","gen_temp		 = 300");	
        fprintf(fp2,"%s\n","constraints         = none");
	fclose(fp2);	
  	 
    system("~/GRO453/bin/grompp -f minim.mdp -c start.gro -p start.top -o run.tpr");

    system("~/GRO453/bin/mdrun -s run.tpr -c start.gro ");

    system("cp ./start.gro ./restart.gro");
    system("cp ./start.top ./restart.top");
   
   return 1;
}

 int solvent_config_separator(int *num_na,int *num_cl,int imp)
 
 {
      
  int i,a,m,a1,k,o,u;  
  
  char char1[100];
      
  char res_type[5],atom_type[2];
  
  int atom_num;
  
  float atom_x,atom_y,atom_z;

  int res_number;  
  
  int natoms,atom_tot;
  
  FILE *fp;
  FILE *fp2;
  
  float box_x,box_y,box_z;
  
  char *chari;
  
  int count_na,count_cl;


  if(imp == 0) {
    
    
      fp = fopen("minimized3.gro","r");
      
      a = 0;
      
      natoms = 0;
      
      count_na = 0;
      
      count_cl = 0; 
      
      while(fgets(char1,sizeof(char1),fp) != NULL){
	
	            a ++;
	
       //	            printf("%s\n",char1);
		    
		if(a ==2) {
		  
		sscanf(char1,"%d",&atom_tot);
		
                };
	
		if(a > 2 && a <= atom_tot+2) {    
	            sscanf(char1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f",&res_number,res_type,atom_type,&atom_num,&atom_x,&atom_y,&atom_z);

	//	    printf("%s\n",atom_type);
	
                    if(strcmp(res_type,"SOL")!=0){
	    
	            natoms++;
	    
	            };

	
                   if(strcmp(atom_type,"NA")==0){
	    
	            count_na++;
	     
	            };	            

                   if(strcmp(atom_type,"CL")==0){
	    
	            count_cl++;
	    
	           };	
		   
		};	
		
		if(a == atom_tot+3) sscanf(char1,"%f%f%f",&box_x,&box_y,&box_z);
	
      };
      
      fclose(fp); 
      
   //   printf("%d\n",count_na);
   //   printf("%d\n",count_cl);
      
   //   printf("%d\t%s\n",natoms,"natoms");
    
  };
  
      *num_na = count_na;
      *num_cl = count_cl;

      
      a = 0;
    
      fp = fopen("minimized3.gro","r");      
      fp2 = fopen("copy.gro","w");
      
      fprintf(fp2,"%s\n","dried");
      fprintf(fp2,"%d\n",natoms);      
      
      while(fgets(char1,sizeof(char1),fp) != NULL){
	
	            a ++;

	
		if(a > 2 && a <= natoms+2) {  
		  
	            sscanf(char1,"%d%s%s%d%f%f%f",&res_number,res_type,atom_type,&atom_num,&atom_x,&atom_y,&atom_z);
                    fprintf(fp2,"%5d%5s%5s%5d%7.3f%7.3f%7.3f\n",res_number,res_type,atom_type,atom_num,atom_x,atom_y,atom_z);
	//	    printf("%s\n",atom_type);
	
		};
		
		if(a == natoms+3) fprintf(fp2,"%f\t%f\t%f\n",box_x,box_y,box_z);		
	
      };
      
      fclose(fp);       
      fclose(fp2);

/*
      
      a = 0;
      
      fp2 = fopen("minimized3.gro","w");      
      fp =  fopen("copy.gro","r");
      
      fprintf(fp2,"%s\n","dried");
      fprintf(fp2,"%d\n",natoms);      
      
      while(fgets(char1,sizeof(char1),fp) != NULL){
	
	            a ++;

	
		if(a > 2 && a <= natoms+2) {  
		  
	            sscanf(char1,"%d%s%s%d%f%f%f",&res_number,res_type,atom_type,&atom_num,&atom_x,&atom_y,&atom_z);
                    fprintf(fp2,"%5d%5s%5s%5d%7.3f%7.3f%7.3f\n",res_number,res_type,atom_type,atom_num,atom_x,atom_y,atom_z);
	//	    printf("%s\n",atom_type);
	
		};
		
		if(a == natoms+3) fprintf(fp2,"%f\t%f\t%f\n",box_x,box_y,box_z);		
	
      };
      
      fclose(fp);       
      fclose(fp2);
*/
    fp2 = fopen("minim.mdp","w");
	
	fprintf(fp2,"%s\n","title                    = implicit minim.");
	fprintf(fp2,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp2,"%s\n","include 		 = -I../top");
        fprintf(fp2,"%s\n","define               = -DPOSRE");
        fprintf(fp2,"%s\n","integrator		 = l-bfgs");
	fprintf(fp2,"%s\n","emtol                = 200");
	fprintf(fp2,"%s\n","emstep               = 0.02");
	fprintf(fp2,"%s\n","nsteps               = 200");
	fprintf(fp2,"%s\n","dt			 = 0.001");
	fprintf(fp2,"%s\n","comm_mode            =  linear ");
	fprintf(fp2,"%s\n","nstxout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstvout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp2,"%s\n","nstenergy		 = 1500 ");
	fprintf(fp2,"%s\n","nstxtcout		 = 1500 ");	
	fprintf(fp2,"%s\n","xtc_grps		 = protein");
	fprintf(fp2,"%s\n","energygrps		 = protein ");	
	fprintf(fp2,"%s\n","nstlist 		 = 10");	
	fprintf(fp2,"%s\n","ns_type 		 = grid");		
	fprintf(fp2,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp2,"%s\n","coulombtype		 = pme");		
	fprintf(fp2,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp2,"%s\n","vdwtype           = shift ");	
	fprintf(fp2,"%s\n","rvdw			 = 1.0");	
	fprintf(fp2,"%s\n","pbc                  = xyz");	
	fprintf(fp2,"%s\n","tcoupl  		 = no");	
	fprintf(fp2,"%s\n","tc-grps 		 = system");	
	fprintf(fp2,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp2,"%s\n","ref_t		 = 300");	
	fprintf(fp2,"%s\n","Pcoupl  		 = no");	
	fprintf(fp2,"%s\n","gen_vel 		 = no");	
	fprintf(fp2,"%s\n","gen_temp		 = 300");	
        fprintf(fp2,"%s\n","constraints         = none");
	fclose(fp2);
      
    system("~/GRO453/bin/grompp -f minim.mdp -c minimized3.gro -p start.top -o run.tpr");

    system("~/GRO453/bin/mdrun -s run.tpr -c minimized3.gro ");
    
  
}

int dihedralbreak(int argc,char *argv[],int scansteps,int aminoresnum,float rot_change,float energy_scan2[1000][7][50],int num_list1,float sim_temp,float OMEGA_B[7],int event_or_scan,int num_event)
{

  int i,a,m,a1,k,o,u;  
  
  float OMEGA_B_VAL, OMEGA_ZERO_VAL;
  
  FILE *fp;
  
  float energy_val;
  
  int int_ex;
 
  
  for(k=1;k<=scansteps;k++){

       
    i = dih_processor(k,aminoresnum,rot_change,scansteps,num_event);    
         

    
    printf("%s\n","barrier reached");
    
    i = energy_shift(argc,argv,&energy_val,sim_temp,int_ex);


    if(int_ex == 1) return 0;

    
    energy_scan2[num_list1][num_event][k] = energy_val;
    
    printf("%f%s\n",energy_scan2[num_list1][num_event][k],"energy");

    
  };

  if(event_or_scan == 0){
	
    i = kramer_search(argc,argv,&OMEGA_B_VAL,sim_temp,aminoresnum);

  };
    
    OMEGA_B[num_event] = OMEGA_B_VAL;
    
 //   OMEGA_ZERO[num_event] = OMEGA_ZERO_VAL;
    
    
}

int dih_processor(int k,int aminoresnum,float rot_change,int scansteps,int num_event){
  
  int i,a,m,a1,o,u;  
  
  char char1[50000][100];
  
  char charA[50000][30];  
      
  char res_type[50000][5],atom_type[50000][5];
  
  int atom_num[50000];
  
  float atom_x[50000],atom_y[50000],atom_z[50000];
  
  float atom_x2[50000],atom_y2[50000],atom_z2[50000];  

  int res_number[50000];
  
  int res_test;
  
  int count1;
  
  float box_x,box_y,box_z;
  
  float box_tester;
  
  float energy_val;
  
  char chari;
  
  float atom_1x,atom_1y,atom_1z;
  
  float atom_2x,atom_2y,atom_2z;
  
  float atom_3x,atom_3y,atom_3z;
  
  float diff_x1,diff_y1,diff_z1;
  
  float n1,n2,n3;
  
  float OMEGA_B_VAL, OMEGA_ZERO_VAL;
  
  FILE *fp;
  
  float PI = 3.141592654;
 
  
    if(k==1){
    
      fp = fopen("minimized3.gro","r");
      
      a = 0;
      
      while(fgets(char1[a],sizeof(char1[a]),fp) != NULL){
	
	a ++;
	
//	printf("%d\n",a);
	
      };
      
      u = 0;
            
      for(i=2;i<=a;i++){
    
	
        if(i<=a-2)  sscanf(char1[i],"%20[ abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890]%f%f%f",&charA[i],&atom_x[i],&atom_y[i],&atom_z[i]);

	   
        if(i < 10000) {

                      sscanf(charA[i],"%d%s%s%d",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i]);
		    
                      };	
	
	if(i==a-1)  {

	            sscanf(char1[i],"%f%f%f%s[^\n]",&box_x,&box_y,&box_z,&chari);
//	            printf("%f%f%f\n",box_x,box_y,box_z);
                    };
	
// 	if(i<=a-2) printf("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",res_number[i],&res_type[i],&atom_type[i],atom_num[i],atom_x[i],atom_y[i],atom_z[i]); 
//      if(i==a-1) printf("%s",char1[i]);  
	
      };   
      
  //     printf("%s",char1[a-3]);
      
    //  sscanf(char1[a-1],"%f%f%f",&box_x,&box_y,&box_z);
      
      fclose(fp);
      
    //  box_tester = box_x[0];
      
   //   printf("%f%f%f\n",box_x[1],box_y[1],box_z[1]);
   //    printf("%d%d%d\n",aminoresnum,scansteps,a);

	fp = fopen("1.ndx","w");
	
	fprintf(fp,"%s\n","[ 1 ]");

        for(m=2;m<=a-2;m++){
	    
	   
	//    printf("%s\n",&atom_type[m]);
	  if(aminoresnum == res_number[m]) {
	    
	        if(strcmp(atom_type[m],"CA")==0)
	    
	           fprintf(fp,"%d\n",atom_num[m]);
	
		if(strcmp(atom_type[m],"N")==0)
		  
		   fprintf(fp,"%d\n",atom_num[m]);
		
		if(strcmp(atom_type[m],"C")==0)
		  
		   fprintf(fp,"%d\n",atom_num[m]);
		
		if(strcmp(atom_type[m],"O")==0)
		  
		   fprintf(fp,"%d\n",atom_num[m]);
	           
	//    if(x[i] = res_number[m] && sscanf(atom_type[m],"%3s%2s",&s1,&s2)==1) printf("%s",&atom_type[m]) ;

	      
	    };
	    
	 };
	  
	fclose(fp);	    

	fp = fopen("2.ndx","w");
	
	fprintf(fp,"%s\n","[ group ]");
	
	for(m=2;m<=a-2;m++){
	  
	  if(res_number[m] == aminoresnum) {
	    
	    fprintf(fp,"%d\n",atom_num[m]);
	    
	    count1++;
	    
	  };
	  
	};
	
	fprintf(fp,"%s\n","[ protein ]");
	
	for(m=2;m<=a-2;m++){
	  
	    fprintf(fp,"%d\n",atom_num[m]);
	};
	
	fprintf(fp,"%s\n","[ system ]");
	
	for(m=2;m<=a-2;m++){
	  
	    fprintf(fp,"%d\n",atom_num[m]);
	};	
	
	fclose(fp);
	
    };
    
    if(k >= 2){
     
      fp = fopen("minimized.gro","r");
      
      a = 0;
      
      while(fgets(char1[a],sizeof(char1[a]),fp) != NULL){
	
	a ++;
	
      };
      
      u = 0;
            
      for(i=2;i<=a;i++){
    
	
        if(i<=a-2)  sscanf(char1[i],"%20[ abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890]%f%f%f",&charA[i],&atom_x[i],&atom_y[i],&atom_z[i]);

	   
        if(i < 10000) {

                      sscanf(charA[i],"%d%s%s%d",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i]);
		    
                      };	
	
	if(i==a-1)  {

	            sscanf(char1[i],"%f%f%f%s[^\n]",&box_x,&box_y,&box_z,&chari);
//	            printf("%f%f%f\n",box_x,box_y,box_z);
                    };
	
// 	if(i<=a-2) printf("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",res_number[i],&res_type[i],&atom_type[i],atom_num[i],atom_x[i],atom_y[i],atom_z[i]); 
//      if(i==a-1) printf("%s",char1[i]);  
	
      }; 
      
    //  sscanf(char1[a-1],"%f%f%f",&box_x,&box_y,&box_z);
      
      fclose(fp);
      
      
    };
    
    count1 = 1;
    
    for(i=2;i<=a;i++){
      
      if(res_number[i] == aminoresnum && count1 == 1 ) {
	
	 atom_1x = atom_x[i+1];
	 atom_1y = atom_y[i+1];
	 atom_1z = atom_z[i+1];
	 
	 atom_2x = atom_x[i+2];
	 atom_2y = atom_y[i+2];
	 atom_2z = atom_z[i+2];
	 
	 atom_3x = atom_x[i+3];
	 atom_3y = atom_y[i+3];
	 atom_3z = atom_z[i+3];	 
	 
	 count1 ++;
	 
      };
      
    };
    
    diff_x1 = atom_2x - atom_3x;
    diff_y1 = atom_2y - atom_3y;
    diff_z1 = atom_2z - atom_3z;
    
    n1 = diff_x1 / sqrt(pow(diff_x1,2)+pow(diff_y1,2)+pow(diff_z1,2));
    n2 = diff_y1 / sqrt(pow(diff_x1,2)+pow(diff_y1,2)+pow(diff_z1,2));
    n3 = diff_z1 / sqrt(pow(diff_x1,2)+pow(diff_y1,2)+pow(diff_z1,2));
    
    printf("%f%f%f\n",n1,n2,n3);
    
    for(i=2;i<=a;i++){
      
      if(res_number[i] == aminoresnum) {
	
          atom_x[i] = atom_x[i] - atom_3x;
          atom_y[i] = atom_y[i] - atom_3y;
          atom_z[i] = atom_z[i] - atom_3z;	  
	
      };
      
    };

    for(i=2;i<=a;i++){
      
      if(res_number[i] == aminoresnum) {      
      
	 atom_x2[i] = atom_x[i]*(cos(rot_change*PI/180) + pow(n1,2)*(1 - cos(rot_change*PI/180))) + atom_y[i]*((n2*n1)*(1 - cos(rot_change*PI/180))
	   	- n3 * sin(rot_change*PI/180)) + atom_z[i]*((n3*n1)*(1 - cos(rot_change*PI/180))+ n2*sin(rot_change*PI/180));

	 atom_y2[i] = atom_x[i]*(n1*n2*(1 - cos(rot_change*PI/180))+n3*(sin(rot_change*PI/180))) + atom_y[i]*( cos(rot_change*PI/180) + pow(n2,2)*	
	        (1 - cos(rot_change*PI/180))) + atom_z[i]*(n3*n1*(1 - cos(rot_change*PI/180))-n1*sin(rot_change*PI/180));
		
	 atom_z2[i] = atom_x[i]*(n1*n3*(1 - cos(rot_change*PI/180))-n2*sin(rot_change*PI/180)) + atom_y[i]*( n2*n3*	(1 - cos(rot_change*PI/180)) +
	        n1*sin(rot_change*PI/180)) + atom_z[i]*( cos(rot_change*PI/180) + pow(n3,2) *(1 - cos(rot_change*PI/180)));
      };
      
    };  

    for(i=2;i<=a;i++){
      
      if(res_number[i] == aminoresnum) {
	
          atom_x2[i] = atom_x2[i] + atom_3x;
          atom_y2[i] = atom_y2[i] + atom_3y;
          atom_z2[i] = atom_z2[i] + atom_3z;	  
	
      };
      
    };    

    fp = fopen("minimized.gro","w");
    
    fprintf(fp,"%s\n","intermediate-struct");
    fprintf(fp,"%d\n", a-3);
    
    for(i=2;i<=a-2;i++){
      
      if(res_number[i] == aminoresnum) {
	
	 atom_x[i] = atom_x2[i];
	 atom_y[i] = atom_y2[i];
	 atom_z[i] = atom_z2[i];
	 
      };
      
      fprintf(fp,"%s%8.3f%8.3f%8.3f\n",charA[i],atom_x[i],atom_y[i],atom_z[i]);
      
    };
    
    fprintf(fp,"%s",char1[a-1]);

    fclose(fp);
    
    printf("%s\n","before energy shift");  
  
    if(k == scansteps){
    
    fp = fopen("1.ndx","w");
	
    fprintf(fp,"%s\n","[ group ]");
	
	for(m=2;m<=a-2;m++){
	  
//	  if(res_number[m] == aminoresnum) {
	    
	    fprintf(fp,"%d\n",atom_num[m]);
	    
	    count1++;
	    
//	  };
	  
	};
	
	fclose(fp);  
	
    };
    
    return 1;
  
}

int kramer_search(int argc,char *argv[],float *OMEGA_B_VAL,float sim_temp,int aminoresnum)
{
    int i,k,j;

    char *grofile,*topfile,*ndxfile;
	 
    char *trajfile,*final_file,*ener_file,*mdpfile;
	
    int mdtime = 10000,num_rls;
    
    int select_term;
    
    float freq_event;
    
    float md_time,stddev_en;
    
    md_time = (float)mdtime;
    
     float rls[10000];

     char *xvgfile;
     int sizeof_x;
     int x[100000];
     float y[100000], sizeof_x2, aver, aver_square;

     FILE *fp2;

    
        FILE *fp = fopen("kramer.mdp","w");
	
	fprintf(fp,"%s\n","title                    = implicit prod.");
	fprintf(fp,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp,"%s\n","include 		 = -I../top");
        fprintf(fp,"%s\n","integrator		 = md");
	fprintf(fp,"%s\n","dt			 = 0.001");
	fprintf(fp,"%s %d\n","nsteps  		 = ",mdtime);
	fprintf(fp,"%s\n","comm_mode         =  linear ");
	fprintf(fp,"%s\n","nstxout 		 = 1500");
	fprintf(fp,"%s\n","nstvout 		 = 1500");
	fprintf(fp,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp,"%s\n","nstenergy		 = 10 ");
	fprintf(fp,"%s\n","nstxtcout		 = 1");	
	fprintf(fp,"%s\n","xtc_grps		 = protein");
	fprintf(fp,"%s\n","energygrps		 = protein ");	
	fprintf(fp,"%s\n","nstlist 		 = 10");	
	fprintf(fp,"%s\n","ns_type 		 = grid");		
	fprintf(fp,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp,"%s\n","coulombtype		 = pme ");		
	fprintf(fp,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp,"%s\n","vdwtype           = shift ");	
	fprintf(fp,"%s\n","rvdw			 = 1.0");	
	fprintf(fp,"%s\n","pbc                      = xyz");	
	fprintf(fp,"%s\n","tcoupl  		 = v-rescale");	
	fprintf(fp,"%s\n","tc-grps 		 = system");	
	fprintf(fp,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp,"%s %f \n","ref_t		 = ",sim_temp);	
	fprintf(fp,"%s\n","Pcoupl  		 = berendsen");
    fprintf(fp,"%s\n","tau_p            =  1.0");
    fprintf(fp,"%s\n","compressibility  =  4E-5");
    fprintf(fp,"%s\n","ref_p            =  1.0");	
	fprintf(fp,"%s\n","gen_vel 		 = yes");	
	fprintf(fp,"%s %f \n","gen_temp		 = ",sim_temp);	
        fprintf(fp,"%s\n","constraints         = none");
//	fprintf(fp,"%s\n","freezegrps               = group");
//	fprintf(fp,"%s\n","freezedim                = Y Y Y");
	fclose(fp);

//	for(i=1;i<=2;i++){
	  
    fp2 = fopen("minim.mdp","w");
	
	fprintf(fp2,"%s\n","title                    = implicit minim.");
	fprintf(fp2,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp2,"%s\n","include 		 = -I../top");
        fprintf(fp2,"%s\n","integrator		 = l-bfgs");
	fprintf(fp2,"%s\n","emtol                = 5000");
	fprintf(fp2,"%s\n","emstep               = 0.02");
	fprintf(fp2,"%s\n","nsteps               = 200");
	fprintf(fp2,"%s\n","dt			 = 0.001");
	fprintf(fp2,"%s\n","comm_mode            =  linear ");
	fprintf(fp2,"%s\n","nstxout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstvout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp2,"%s\n","nstenergy		 = 1500 ");
	fprintf(fp2,"%s\n","nstxtcout		 = 1500 ");	
	fprintf(fp2,"%s\n","xtc_grps		 = protein");
	fprintf(fp2,"%s\n","energygrps		 = protein ");	
	fprintf(fp2,"%s\n","nstlist 		 = 10");	
	fprintf(fp2,"%s\n","ns_type 		 = grid");		
	fprintf(fp2,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp2,"%s\n","coulombtype		 = pme");		
	fprintf(fp2,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp2,"%s\n","vdwtype           = shift ");	
	fprintf(fp2,"%s\n","rvdw			 = 1.0");	
	fprintf(fp2,"%s\n","pbc                  = xyz");	
	fprintf(fp2,"%s\n","tcoupl  		 = no");	
	fprintf(fp2,"%s\n","tc-grps 		 = system");	
	fprintf(fp2,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp2,"%s\n","ref_t		 = 300");	
	fprintf(fp2,"%s\n","Pcoupl  		 = no");	
	fprintf(fp2,"%s\n","gen_vel 		 = no");	
	fprintf(fp2,"%s\n","gen_temp		 = 300");	
        fprintf(fp2,"%s\n","constraints         = none");
	fclose(fp2);
      
    system("~/GRO453/bin/grompp -f minim.mdp -c minimized.gro -p start.top -o run.tpr");

    system("~/GRO453/bin/mdrun -s run.tpr -c kramer.gro ");
	
        system("cp LOV2.top kramer.top");
	 
 //       if(i == 1) system("~/gromacs-4.5.5/src/tools/genbox -cp minimized3.gro -cs spc216.gro -p kramer.top -o kramer.gro");
 //    system("~/gromacs-4.5.5/src/tools/genbox -cp minimized.gro -cs spc216.gro -p kramer.top -o kramer.gro");   
 
    fp2 = fopen("minim.mdp","w");
	
	fprintf(fp2,"%s\n","title                    = implicit minim.");
	fprintf(fp2,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp2,"%s\n","include 		 = -I../top");
        fprintf(fp2,"%s\n","integrator		 = cg");
	fprintf(fp2,"%s\n","emtol                = 2000");
	fprintf(fp2,"%s\n","emstep               = 0.02");
	fprintf(fp2,"%s\n","nsteps               = 200");
	fprintf(fp2,"%s\n","dt			 = 0.001");
	fprintf(fp2,"%s\n","comm_mode            =  linear ");
	fprintf(fp2,"%s\n","nstxout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstvout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp2,"%s\n","nstenergy		 = 1500 ");
	fprintf(fp2,"%s\n","nstxtcout		 = 1500 ");	
	fprintf(fp2,"%s\n","xtc_grps		 = protein");
	fprintf(fp2,"%s\n","energygrps		 = protein ");	
	fprintf(fp2,"%s\n","nstlist 		 = 10");	
	fprintf(fp2,"%s\n","ns_type 		 = grid");		
	fprintf(fp2,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp2,"%s\n","coulombtype		 = pme");		
	fprintf(fp2,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp2,"%s\n","vdwtype           = shift ");	
	fprintf(fp2,"%s\n","rvdw			 = 1.0");	
	fprintf(fp2,"%s\n","pbc                  = xyz");	
	fprintf(fp2,"%s\n","tcoupl  		 = no");	
	fprintf(fp2,"%s\n","tc-grps 		 = system");	
	fprintf(fp2,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp2,"%s\n","ref_t		 = 300");	
	fprintf(fp2,"%s\n","Pcoupl  		 = no");	
	fprintf(fp2,"%s\n","gen_vel 		 = no");	
	fprintf(fp2,"%s\n","gen_temp		 = 300");	
        fprintf(fp2,"%s\n","constraints         = none");
	fclose(fp2);	
  	 
 //   system("~/gromacs-4.5.5/src/kernel/grompp -f minim.mdp -c kramer.gro -p kramer.top -o run.tpr");

//   system("mpirun -np 8 ~/gromacs-4.5.5/src/kernel/mdrun -s run.tpr -c kramer.gro");
     
//        grompp_float(argc,argv);
        system("~/GRO453/bin/grompp -f kramer.mdp -p start.top -c kramer.gro -maxwarn 3 -o run.tpr ");	
	
	system("mpirun -np 6 ~/gromacs-4.5.5/src/kernel/mdrun -v -s run.tpr -c kramer.gro -o trj.trr -e md_ener.edr -pd");

	
//	gmx_rms2(argc,argv,&num_rls,rls);

        select_term = 13;
	
	
    //    gmx_energy2(argc,argv,select_term,&stddev_en);
	
    fp = fopen("en.sh","w");
    fprintf(fp,"%s\n","#!/bin/bash");
    fprintf(fp,"%s\n","~/GRO453/bin/g_energy -f md_ener.edr -o t.xvg << EOF");
    fprintf(fp,"%s\n","13");
    fprintf(fp,"%s\n","0");
    fprintf(fp,"%s\n","EOF");
    fclose(fp);

    system("chmod 744 ./en.sh");

    system("./en.sh");


    xvgfile = "t.xvg";    
 
    read_xvg_files(xvgfile,&sizeof_x,x,y);

   /* for(k=1;k<=sizeof_x;k++){

          printf("%g\n",y[k]);

    }; */

    aver = 0;
    aver_square = 0;

    for(k=1;k<=sizeof_x;k++){

                 aver = aver + y[k];
          aver_square = aver_square + pow(y[k],2);

    };

    if(sizeof_x == 0) {

       sizeof_x = 1;
       aver     = 1;
    aver_square = 1;

    };

    sizeof_x2 = (float)sizeof_x;

    aver = aver / sizeof_x2;

    aver_square = aver_square / sizeof_x2;

    stddev_en = sqrt(aver_square - pow(aver,2));

	printf("%g\t%s\n",stddev_en,"STDDEVEBN");
	
	
  
	  *OMEGA_B_VAL = stddev_en/md_time;

	if(isnanf(*OMEGA_B_VAL) == 1) *OMEGA_B_VAL = 1;	  	
//	};
  
	printf("%e\n",*OMEGA_B_VAL);
 	
	
}

int kramer_zero(int argc,char *argv[],float *OMEGA_B_VAL,float sim_temp,int aminoresnum)
{
    int i,k,j;

    char *grofile,*topfile,*ndxfile;
	 
    char *trajfile,*final_file,*ener_file,*mdpfile;
	
    int mdtime = 10000,num_rls;
    
    int select_term;
    
    float freq_event;
    
    float md_time,stddev_en;
    
    md_time = (float)mdtime;
    
     float rls[10000];

     char *xvgfile;
     int sizeof_x;
     int x[100000];
     float y[100000], sizeof_x2, aver, aver_square;

     FILE *fp2;

     FILE *fp = fopen("kramer.mdp","w");
	
	fprintf(fp,"%s\n","title                    = implicit prod.");
	fprintf(fp,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp,"%s\n","include 		 = -I../top");
    fprintf(fp,"%s\n","integrator		 = md");
	fprintf(fp,"%s\n","dt			 = 0.001");
	fprintf(fp,"%s %d\n","nsteps  		 = ",mdtime);
	fprintf(fp,"%s\n","comm_mode         =  linear ");
	fprintf(fp,"%s\n","nstxout 		 = 1500");
	fprintf(fp,"%s\n","nstvout 		 = 1500");
	fprintf(fp,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp,"%s\n","nstenergy		 = 10 ");
	fprintf(fp,"%s\n","nstxtcout		 = 1");	
	fprintf(fp,"%s\n","xtc_grps		 = protein");
	fprintf(fp,"%s\n","energygrps		 = protein ");	
	fprintf(fp,"%s\n","nstlist 		 = 10");	
	fprintf(fp,"%s\n","ns_type 		 = grid");		
	fprintf(fp,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp,"%s\n","coulombtype		 = pme ");		
	fprintf(fp,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp,"%s\n","vdwtype           = shift ");	
	fprintf(fp,"%s\n","rvdw			 = 1.0");	
	fprintf(fp,"%s\n","pbc                      = xyz");	
	fprintf(fp,"%s\n","tcoupl  		 = v-rescale");	
	fprintf(fp,"%s\n","tc-grps 		 = system");	
	fprintf(fp,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp,"%s %f \n","ref_t		 = ",sim_temp);	
	fprintf(fp,"%s\n","Pcoupl  		 = berendsen");
        fprintf(fp,"%s\n","tau_p            =  1.0");
        fprintf(fp,"%s\n","compressibility  =  4E-5");
        fprintf(fp,"%s\n","ref_p            =  1.0");	
	fprintf(fp,"%s\n","gen_vel 		 = yes");	
	fprintf(fp,"%s %f \n","gen_temp		 = ",sim_temp);	
        fprintf(fp,"%s\n","constraints         = none");   
//	fprintf(fp,"%s\n","freezegrps               = group");
//      fprintf(fp,"%s\n","freezedim                = Y Y Y");		
	fclose(fp);

//	for(i=1;i<=2;i++){
    fp2 = fopen("minim.mdp","w");
	
	fprintf(fp2,"%s\n","title                    = implicit minim.");
	fprintf(fp2,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp2,"%s\n","include 		 = -I../top");
        fprintf(fp2,"%s\n","integrator		 = cg");
	fprintf(fp2,"%s\n","emtol                = 5000");
	fprintf(fp2,"%s\n","emstep               = 0.02");
	fprintf(fp2,"%s\n","nsteps               = 100");
	fprintf(fp2,"%s\n","dt			 = 0.001");
	fprintf(fp2,"%s\n","comm_mode            =  linear ");
	fprintf(fp2,"%s\n","nstxout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstvout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp2,"%s\n","nstenergy		 = 1500 ");
	fprintf(fp2,"%s\n","nstxtcout		 = 1500 ");	
	fprintf(fp2,"%s\n","xtc_grps		 = protein");
	fprintf(fp2,"%s\n","energygrps		 = protein ");	
	fprintf(fp2,"%s\n","nstlist 		 = 10");	
	fprintf(fp2,"%s\n","ns_type 		 = grid");		
	fprintf(fp2,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp2,"%s\n","coulombtype		 = pme");		
	fprintf(fp2,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp2,"%s\n","vdwtype           = shift ");	
	fprintf(fp2,"%s\n","rvdw			 = 1.0");	
	fprintf(fp2,"%s\n","pbc                  = xyz");	
	fprintf(fp2,"%s\n","tcoupl  		 = no");	
	fprintf(fp2,"%s\n","tc-grps 		 = system");	
	fprintf(fp2,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp2,"%s\n","ref_t		 = 300");	
	fprintf(fp2,"%s\n","Pcoupl  		 = no");	
	fprintf(fp2,"%s\n","gen_vel 		 = no");	
	fprintf(fp2,"%s\n","gen_temp		 = 300");	
        fprintf(fp2,"%s\n","constraints         = none");
	fclose(fp2);
      
    system("~/GRO453/bin/grompp -f minim.mdp -c minimized3.gro -p start.top -o run.tpr");

    system("~/GRO453/bin/mdrun -s run.tpr -c minimized3.gro ");	  
	
        system("cp LOV2.top kramer.top");
	system("cp minimized3.gro kramer.gro");
	 
 //       if(i == 1) system("~/gromacs-4.5.5/src/tools/genbox -cp minimized3.gro -cs spc216.gro -p kramer.top -o kramer.gro");
 //       system("~/gromacs-4.5.5/src/tools/genbox -cp minimized3.gro -cs spc216.gro -p kramer.top -o kramer.gro");    

    fp2 = fopen("minim.mdp","w");
	
	fprintf(fp2,"%s\n","title                    = implicit minim.");
	fprintf(fp2,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp2,"%s\n","include 		 = -I../top");
        fprintf(fp2,"%s\n","integrator		 = cg");
	fprintf(fp2,"%s\n","emtol                = 2000");
	fprintf(fp2,"%s\n","emstep               = 0.02");
	fprintf(fp2,"%s\n","nsteps               = 200");
	fprintf(fp2,"%s\n","dt			 = 0.001");
	fprintf(fp2,"%s\n","comm_mode            =  linear ");
	fprintf(fp2,"%s\n","nstxout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstvout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp2,"%s\n","nstenergy		 = 1500 ");
	fprintf(fp2,"%s\n","nstxtcout		 = 1500 ");	
	fprintf(fp2,"%s\n","xtc_grps		 = protein");
	fprintf(fp2,"%s\n","energygrps		 = protein ");	
	fprintf(fp2,"%s\n","nstlist 		 = 10");	
	fprintf(fp2,"%s\n","ns_type 		 = grid");		
	fprintf(fp2,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp2,"%s\n","coulombtype		 = pme");		
	fprintf(fp2,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp2,"%s\n","vdwtype           = shift ");	
	fprintf(fp2,"%s\n","rvdw			 = 1.0");	
	fprintf(fp2,"%s\n","pbc                  = xyz");	
	fprintf(fp2,"%s\n","tcoupl  		 = no");	
	fprintf(fp2,"%s\n","tc-grps 		 = system");	
	fprintf(fp2,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp2,"%s\n","ref_t		 = 300");	
	fprintf(fp2,"%s\n","Pcoupl  		 = no");	
	fprintf(fp2,"%s\n","gen_vel 		 = no");	
	fprintf(fp2,"%s\n","gen_temp		 = 300");	
        fprintf(fp2,"%s\n","constraints         = none");
	fclose(fp2);
	  	 
//    system("~/gromacs-4.5.5/src/kernel/grompp -f minim.mdp -c kramer.gro -p kramer.top -o run.tpr");

//    system("mpirun -np 8 ~/gromacs-4.5.5/src/kernel/mdrun -s run.tpr -c kramer.gro");     
//        grompp_float(argc,argv);
        system("~/GRO453/bin/grompp -f kramer.mdp -p start.top -c kramer.gro -maxwarn 3 -o run.tpr ");	
	
	system("mpirun -np 6 ~/gromacs-4.5.5/src/kernel/mdrun -v -s run.tpr -c kramer.gro -o trj.trr -e md_ener.edr -pd ");

	
//	gmx_rms2(argc,argv,&num_rls,rls);

        select_term = 13;
	
	
    //    gmx_energy2(argc,argv,select_term,&stddev_en);

    fp = fopen("en.sh","w");
    fprintf(fp,"%s\n","#!/bin/bash");
    fprintf(fp,"%s\n","~/GRO453/bin/g_energy -f md_ener.edr -o t.xvg << EOF");
    fprintf(fp,"%s\n","13");
    fprintf(fp,"%s\n","0");
    fprintf(fp,"%s\n","EOF");
    fclose(fp);

    system("chmod 744 ./en.sh");

    system("./en.sh");

	
    xvgfile = "t.xvg";    
 
    read_xvg_files(xvgfile,&sizeof_x,x,y);

   /* for(k=1;k<=sizeof_x;k++){

          printf("%g\n",y[k]);

    }; */

    aver = 0;
    aver_square = 0;

    for(k=1;k<=sizeof_x;k++){

                 aver = aver + y[k];
          aver_square = aver_square + pow(y[k],2);

    };

    if(sizeof_x == 0) {

       sizeof_x = 1;
       aver     = 1;
   aver_square  = 1;

    };

    sizeof_x2 = (float)sizeof_x;

    aver = aver / sizeof_x2;

    aver_square = aver_square / sizeof_x2;

    stddev_en = sqrt(aver_square - pow(aver,2));

	printf("%g\t%s\n",stddev_en,"STDDEVEBN");
	
	
  
	  *OMEGA_B_VAL = stddev_en/md_time;

	if(isnanf(*OMEGA_B_VAL) == 1) *OMEGA_B_VAL = 1;
	  
//	};
  
	printf("%e\n",*OMEGA_B_VAL);
 	
	
}


int noe_reader(char *xvgfile,int *sizeof_x,int x[100000],float y[100000])
{
  
  int i,a,u;
  
  int step;
  
  float value;
  
  char char1 [50];
  
  char s[1],c[49];
    
  
       FILE *fp = fopen(xvgfile,"r");
       
       u = 1;
       a = 1;
       	
        while(fgets(char1,sizeof(char1),fp) != 0) 
	{
	    
	     sscanf(char1,"%d%f",&step,&value);
	 
	     x[u] = step;
	     
	     y[u] = value;
	     
	     *sizeof_x = u;
	    
	    printf("%f\t%f\t%s\n",step,value,"value");
        //    printf("%s\n",char1);	 
 
	     u ++;	 
	  
	};
      
       	printf("%d\n",u);
	
//	u = 1;
	
/*	for(i = 0; i <= a; i++)
	{
	  if(sscanf(char1,"#%s",s,c)==0 && sscanf(char1,"@%s",s,c)==0) { 
	    
	     sscanf(char1,"%d%f",&step,&value);
	 
	     x[u] = step;
	     
	     y[u] = value;
	     
	     *sizeof_x = u;
	    
	 //    printf("%d\t%f\n",step[u],value[u]);
	  
	     u ++;
	     
	  };
	  
	    
	 
	};*/  
	
	return 0; 

}

int energy_shift(int argc,char *argv[],float *energy_val,float sim_temp,int int_ex)

{
  int i,k;
  
  FILE *fp;
  
  float ran;
  
  float lambda;
  
  char *grofile, *topfile, *mdpfile, *ndxfile;  

  char *trajfile;
  
  char *final_file;

  char *ener_file;

  char *dhdlfile;
  
  float free_energy_out;
  
  i = random_number(&ran);
  
  ran = ran * 100000;
  
  lambda = 0;

   
  for(k = 1;k <= 3; k++) {  
    

    fp = fopen("minim_fee.mdp","w");

    fprintf(fp,"%s\n","include                  = -I../top");
    fprintf(fp,"%s\n","integrator               = md");
    fprintf(fp,"%s\n","tinit                    = 0");
    fprintf(fp,"%s\n","dt                       = 0.0005");
    fprintf(fp,"%s\n","nsteps                   = 100");
    fprintf(fp,"%s\n","comm-mode                = Linear");
    fprintf(fp,"%s\n","nstxout                  = 1");
    fprintf(fp,"%s\n","nstvout                  = 1");
    fprintf(fp,"%s\n","nstfout                  = 0");
    fprintf(fp,"%s\n","nstlog                   = 1");
    fprintf(fp,"%s\n","nstcalcenergy            = 5");
    fprintf(fp,"%s\n","nstenergy                = 5");
    fprintf(fp,"%s\n","nstxtcout                = 10");
    fprintf(fp,"%s\n","xtc-precision            = 1000");
    fprintf(fp,"%s\n","xtc-grps                 = protein");
    fprintf(fp,"%s\n","energygrps               = protein");
    fprintf(fp,"%s\n","nstlist                  = 1");
    fprintf(fp,"%s\n","ns_type                  = grid");
    fprintf(fp,"%s\n","pbc                      = xyz");
    fprintf(fp,"%s\n","rlist                    = 1.0");
    fprintf(fp,"%s\n","coulombtype              = shift");
    fprintf(fp,"%s\n","rcoulomb-switch          = 0");
    fprintf(fp,"%s\n","rcoulomb                 = 1.0");
    fprintf(fp,"%s\n","vdwtype                  = shift");      
    fprintf(fp,"%s\n","rvdw-switch              = 0");
    fprintf(fp,"%s\n","rvdw                     = 1.0");
    fprintf(fp,"%s\n","table-extension          = 100"); 
    fprintf(fp,"%s\n","pme_order                = 4");
    fprintf(fp,"%s\n","ewald_rtol               = 1e-05");
    fprintf(fp,"%s\n","ewald_geometry           = 3d");
    fprintf(fp,"%s\n","epsilon_surface          = 0");
    fprintf(fp,"%s\n","optimize_fft             = no");
    fprintf(fp,"%s\n","tcoupl                   = v-rescale");
    fprintf(fp,"%s\n","tc-grps                  = system ");
    fprintf(fp,"%s\n","tau-t                    = 1.0 ");
    fprintf(fp,"%s%f\n","ref-t                    = ", sim_temp);
    fprintf(fp,"%s\n","gen-vel                  = yes");
    fprintf(fp,"%s%f\n","gen-temp                 = ",sim_temp);
    fprintf(fp,"%s%d\n","gen-seed                 = ", (int)ran);
    fprintf(fp,"%s\n","constraints              = none");
    fprintf(fp,"%s\n","free_energy              = yes");
    fprintf(fp,"%s%f\n","init_lambda              = ",lambda);
    fprintf(fp,"%s\n","delta_lambda             = 0"); 
    fprintf(fp,"%s\n","foreign_lambda           = 0 0.5 1");
    fprintf(fp,"%s\n","sc-alpha                 = 0.5");
    fprintf(fp,"%s\n","sc-power                 = 1");
    fprintf(fp,"%s\n","sc-sigma                 = 0.3");
    fprintf(fp,"%s\n","couple-moltype           = Protein");
    fprintf(fp,"%s\n","couple-lambda0           = vdw-q");
    fprintf(fp,"%s\n","couple-lambda1           = vdw");
    fprintf(fp,"%s\n","couple-intramol          = yes");
    fprintf(fp,"%s\n","nstdhdl                  = 5");
//    fprintf(fp,"%s\n","freezegrps               = group ");
//    fprintf(fp,"%s\n","freezedim                = Y Y Y ");
//    fprintf(fp,"%s\n","implicit_solvent         = GBSA");
//    fprintf(fp,"%s\n","gb_algorithm             = HCT");
//    fprintf(fp,"%s\n","rgbradii                 = 0.8");

    fclose(fp);

    fp = fopen("minim_dih.mdp","w"); 

    fprintf(fp,"%s\n","cpp             = /lib/cpp");
    fprintf(fp,"%s\n","include         = -I../top  ");
    fprintf(fp,"%s\n","integrator      = l-bfgs");    
    fprintf(fp,"%s\n","emstep          = 0.005");
    fprintf(fp,"%s\n","emtol           = 5000 ");
    fprintf(fp,"%s\n","nsteps          = 100");   
    fprintf(fp,"%s\n","nstenergy       = 1");
    fprintf(fp,"%s\n","nstxtcout       = 10");    
    fprintf(fp,"%s\n","xtc_grps        = Protein");    
    fprintf(fp,"%s\n","energygrps      = Protein");
    fprintf(fp,"%s\n","nstlist         = 5");    
    fprintf(fp,"%s\n","ns_type         = grid");
    fprintf(fp,"%s\n","rlist           = 1.0");   
    fprintf(fp,"%s\n","coulombtype     = pme"); 
    fprintf(fp,"%s\n","rcoulomb        = 1.0");    
    fprintf(fp,"%s\n","vdwtype         = shift"); 
    fprintf(fp,"%s\n","rvdw            = 1.0");    
    fprintf(fp,"%s\n","constraints     = none");   
    fprintf(fp,"%s\n","pbc             = xyz");
//    fprintf(fp,"%s\n","freezegrps      = group ");    
//    fprintf(fp,"%s\n","freezedim       = Y Y Y ");
    

    fclose(fp);
    
    fp = fopen("minim.mdp","w"); 

    fprintf(fp,"%s\n","cpp             = /lib/cpp");
    fprintf(fp,"%s\n","include         = -I../top  ");
    fprintf(fp,"%s\n","integrator      = l-bfgs");    
    fprintf(fp,"%s\n","emstep          = 0.02");
    fprintf(fp,"%s\n","emtol           = 7000 ");
    fprintf(fp,"%s\n","nsteps          = 10");   
    fprintf(fp,"%s\n","nstenergy       = 1");
    fprintf(fp,"%s\n","nstxtcout       = 10");    
    fprintf(fp,"%s\n","xtc_grps        = Protein");    
    fprintf(fp,"%s\n","energygrps      = Protein");
    fprintf(fp,"%s\n","nstlist         = 5");    
    fprintf(fp,"%s\n","ns_type         = simple");
    fprintf(fp,"%s\n","rlist           = 1.0");   
    fprintf(fp,"%s\n","coulombtype     = pme"); 
    fprintf(fp,"%s\n","rcoulomb        = 1.0");    
    fprintf(fp,"%s\n","vdwtype         = shift"); 
    fprintf(fp,"%s\n","rvdw            = 1.0");    
    fprintf(fp,"%s\n","constraints     = none");   
    fprintf(fp,"%s\n","pbc             = xyz");
    
    fclose(fp);    

    lambda = lambda + 0.5;
    
//    system("mv ./2.ndx ./8.ndx");
    
//    system("cat 8.ndx freeze.ndx > 2.ndx");
    
    printf("%s\n","before minimize struct");
    
    i = minimize_struct(argc,argv,int_ex);  
       
  
       system("~/GRO453/bin/grompp -f minim_fee.mdp -p start.top -c minimized.gro -maxwarn 3 -o run.tpr ");    
     
       
 if(k == 1)   system("~/GRO453/bin/mdrun  -v -s run.tpr -c minimized.gro -o trj.trr -e ener.edr -dhdl dhdl1.xvg");   
 if(k == 2)   system("~/GRO453/bin/mdrun  -v -s run.tpr -c minimized.gro -o trj.trr -e ener.edr -dhdl dhdl2.xvg"); 
 if(k == 3)   system("~/GRO453/bin/mdrun  -v -s run.tpr -c minimized.gro -o trj.trr -e ener.edr -dhdl dhdl3.xvg");
 
    system("rm -f ./#*#");

    system("rm -f ./core*");

  }  


   //  i = gmx_bar3(argc,argv,&free_energy_out);
 
    i = barint(&free_energy_out);

    *energy_val = free_energy_out;

}

int minimize_struct(int argc,char *argv[],int int_ex)
{
  
  int i,k;
  
  FILE *fp;
  
  float ran;
  
  float lambda;
  
  char *grofile, *topfile, *mdpfile, *ndxfile;  

  char *trajfile;
  
  char *final_file;

  char *ener_file;

  char *dhdlfile;
  
  float free_energy_out; 

  int sizeof_x;

  float sizeof_xf;

  float x[100000], y[100000];
 
  float average;

  char *xvgfile;


   int_ex = 0;
    
//    i = grompp_float(argc,argv);
    system("~/GRO453/bin/grompp -f minim_dih.mdp -p start.top -c minimized.gro  -maxwarn 3 -o run.tpr");     
	
    system("~/GRO453/bin/mdrun -v -s run.tpr -c minimized.gro -o trj.trr -e ener.edr ");
   
/*    fp = fopen("en.sh","w");
    fprintf(fp,"%s\n","#!/bin/bash");
    fprintf(fp,"%s\n","g_energy -f md_ener.edr -o t.xvg << EOF");
    fprintf(fp,"%s\n","13");
    fprintf(fp,"%s\n","0");
    fprintf(fp,"%s\n","EOF");
    fclose(fp);

    system("chmod 744 ./en.sh");

    system("./en.sh");

 
    xvgfile = "t.xvg";

    read_xvg_files(xvgfile,&sizeof_x,x,y);  
 
    average = 0;

    for(i=1;i<=sizeof_x;i++){

        average = average + y[i];

    }; 

    sizeof_xf = (float)sizeof_x;

    average = average / sizeof_xf;

    if(average > 100000) {

           int_ex = 1;

    }; */

    int_ex = 0;

    printf("%s\n","after minimize mdfloat");  
}


int dih_trans(int argc,char *argv[],float ang_stddev[10000],int *sizeof_stddev)
{
      int i,m,a,a1,k,count1;
  
      char *xvgfile;
        
      int sizeof_x;
      
      int sizeof_x_searchlist;
      
      int natoms;
      
      char char1[50000][100];
      
      char res_type[50000],atom_type[50000][5];
  
      int atom_num[50000];
  
      float atom_x[50000],atom_y[50000],atom_z[50000];

      int res_number[50000]; 
      
      char s1[2];
      
      char s2[3];
      
      char s3;
      
      float angle_stddev;

      float x[100000];

      float y[100000];

      int   x_searchlist[50000];
      
      float y_searchlist[50000];
      
      float sizeof_x2;   

      float aver,aver_square;   
//      float ang_stddev[1000];
      
      FILE *fp;

      
      xvgfile = "searchlist.ndx";    
 
      i = read_xvg_files_int(xvgfile,&sizeof_x_searchlist,x_searchlist,y_searchlist);
     
      
      fp = fopen("minimized3.gro","r");
      
      a = 0;
      
      while(fgets(char1[a],sizeof(char1),fp) != NULL){
	
	a ++;
	
      };
      
      for(i=2;i<=a-2;i++){
    
          sscanf(char1[i],"%d%s%s%d%f%f%f",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i],&atom_x[i],&atom_y[i],&atom_z[i]);
      //    printf("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",res_number[i],&res_type[i],&atom_type[i],atom_num[i],atom_x[i],atom_y[i],atom_z[i]); 
      };      
      
      fclose(fp);
     
 
a1 = 0;      

count1 = 0;

      for(i=1;i<=sizeof_x_searchlist;i++){

	fp = fopen("1.ndx","w");
	
	fprintf(fp,"%s\n","[ 1 ]");
	
	  for(m=2;m<=a-2;m++){
	  
	    if(x_searchlist[i] == res_number[m]) {
	    
	        if(strcmp(atom_type[m],"CA")==0)
	    
	           fprintf(fp,"%d\n",atom_num[m]);
	
		   count1 ++;
		
		if(strcmp(atom_type[m],"N")==0)
		  
		   fprintf(fp,"%d\n",atom_num[m]);
		
		   count1 ++;		
		
		if(strcmp(atom_type[m],"C")==0)
		  
		   fprintf(fp,"%d\n",atom_num[m]);

		   count1 ++;		
		
		if(strcmp(atom_type[m],"O")==0)
		  
		   fprintf(fp,"%d\n",atom_num[m]);

		   count1 ++;		
		
	//    if(x[i] = res_number[m] && sscanf(atom_type[m],"%3s%2s",&s1,&s2)==1) printf("%s",&atom_type[m]) ;

	      
	    };
	    
	  };
	
	if(count1 == 3){
	  
	  fprintf(fp,"%d",atom_num[a-2]);
	  
	};

	
	fclose(fp);

	/*
	argc = 7;
	argv[0] = "g_angle";
	argv[1] = "-f";
	argv[2] = "md_traj.trr";
	argv[3] = "-n";
	argv[4] = "1.ndx";
	argv[5] = "-type";
	argv[6] = "dihedral";
	*/
	  
//	k = gmx_g_angle2(argc,argv,&angle_stddev);

    system("~/GRO453/bin/g_angle -f md_traj.trr -n 1.ndx -type dihedral -ov angle.xvg");

    xvgfile = "angle.xvg";    
 
    read_xvg_files(xvgfile,&sizeof_x,x,y);

    aver = 0;
    aver_square = 0;
    
    if(sizeof_x != 0){

    for(k=1;k<=sizeof_x;k++){

                 aver = aver + y[k];
          aver_square = aver_square + pow(y[k],2);

    };

    sizeof_x2 = (float)sizeof_x;

    aver = aver / sizeof_x2;

    aver_square = aver_square / sizeof_x2;

    angle_stddev = sqrt(aver_square - pow(aver,2));
	
	system("rm -f ./#*");
	
	ang_stddev[i] = angle_stddev;
	
      };
      
    if(sizeof_x == 0) {
      
        ang_stddev[i] = 0.25;
        
      };    
     
      };
      
      *sizeof_stddev = sizeof_x_searchlist;
     
      
    /*  
      for(i=1;i<=sizeof_x_searchlist;i++) {
	
	printf("%f\n",ang_stddev[i]);
	
      }; */
  
}


int noe_sep(int argc,char *argv[],int range_checker_A,int range_checker_B
, int eventmax)
{
  
       int a,i,k,l,p;
  
       int range_1or2;
       
       char *trajfile;
             
       int sizeof_x,sizeof_y;
       
       char *xvgfile;
       
       int x[100000];
       
       float y[100000];
       
       float ry[500000];
       
       float rx[500000];
       
       float maxi_king,mini_king;
       
       int sizerms;
       
       FILE *fp,*fp2;
      
       char res_type[5],atom_type[5];
  
       int atom_num[10000];
  
       float atom_x,atom_y,atom_z;

       int res_number[10000];
       
      float ran, num_class_real_f ;

      int   selected_class;

      int   class_pupils[11][100];

      float pupil_rms[11][100];
      
      float rms_res[11][100];      

      int   class_res [11][100];      

      float test_rank[1000],kf;
      
      int  number_of_pupils[11], empty_class[11];
      
      float class_rank;
      
      int o,n,numat,class_num[1000],number_of_residues,atom_tot;
    
      float numat_f,number_of_residues_f;
 
      float aver_sq, aver;

      float sizerms_f;

  char charA[50000][30];
 
  char char1[50000][100];  
  
      FILE *ndx;
      
      FILE *outkmc;
       
/*    fp = fopen("copy.gro","r");
      fp2 = fopen("prot.ndx","w");
      
      fprintf(fp2,"%s\n","[ 1 ]");
      
      a = 0;
      
      while(fgets(char1,sizeof(char1),fp) != NULL){
	
	a ++;
	
      };
      
      for(i=2;i<=a-2;i++){
    
              fprintf(fp2,"%d\n",i-1);
      };      
      
      fclose(fp); 
      fclose(fp2); */
      
      numat = a-2;
      
      trajfile = "md_traj";
       
       range_1or2 = 1;
       
   /*    argc = 10;
       argv[0] = "g_rmsf";
       argv[1] = "-f";
       argv[2] = "md_traj.trr";
       argv[3] = "-o";
       argv[4] = "rmsf.xvg";
       argv[5] = "-res";
       argv[6] = "-n";
       argv[7] = "prot.ndx";
       argv[8] = "-s";
       argv[9] = "run.tpr"; */
        
   //    i = gmx_rmsf2(argc,argv,range_1or2,trajfile,rmsf2,&sizerms);
 
         fp = fopen("minimized3.gro","r");
      
          a = 0;
      
          while(fgets(char1[a],sizeof(char1[a]),fp) != NULL){
	
        	a ++;
	
          };
            
          for(i=2;i<=a;i++){
    
	
           if(i<=a-2) sscanf(char1[i],"%20[ abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890]%f%f%f",&charA[i],&atom_x,&atom_y,&atom_z);
	     // sscanf(char1[i],"%5d%-5s%5s%5d%8.3f%8.3f%8.3f",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i],&atom_x[i],&atom_y[i],&atom_z[i]);
	   
           if(i < 10000) {

                      sscanf(charA[i],"%d%s%s%d",&res_number[i],&res_type[i],&atom_type[i],&atom_num[i]);
		    
                         };	   
	   
			 
		     if(res_number[i] != res_number[i-1]) number_of_residues = number_of_residues + 1;			 
 	
         };   

         fclose(fp);
	 
      numat  = a-2;
	 
     numat_f = (float) atom_tot;     

number_of_residues_f = (float) number_of_residues;

      eventmax = 5;
 
      for(i=range_checker_A ; i<=range_checker_B ; i++){
	
           ndx = fopen("prot.ndx","w");     

	   n = 0;
	   
	 for(k=1;k<10000;k++){  
	   
           if(res_number[k] == i && n == 0) { 
	     
	     fprintf(ndx,"%s\n","[ 1 ]");
	     
	     n ++;
	  };
		    
	   if(res_number[k] == i) fprintf(ndx,"%d\n", atom_num[k]);  
	   
	 };
	 
	 fclose(ndx);
	 	 
         system("~/GRO453/bin/g_rms -f md_traj.trr -o rms.xvg -n prot.ndx -s minimized3.gro "); 
  
	 system("rm -f ./#*#");
	 
         read_xvg_files("rms.xvg",&sizerms,rx,ry);
	 
	 x[i] = i;
	  
	// y[i] = ry[1] - ry[sizerms];

         aver    = 0;
         aver_sq = 0;

         for(k=1;k<=sizerms;k++){

	    aver = aver + ry[k]; 
         aver_sq = aver_sq + pow(ry[k],2);

         };

         sizerms_f = (float) sizerms;

         aver_sq = aver_sq/sizerms_f;

         aver    = aver / sizerms_f;

         y[i] = sqrt(aver_sq - pow(aver,2));

	// if(y[i] < 0) y[i] = -ry[sizerms] ;//- ry[1];
	 
//	 printf("%f\t%f\n",ry[1],ry[sizerms]);
	 
      };
      
      sizerms = range_checker_B - range_checker_A+1;
      
      /* for(i=1;i<=sizerms;i++){
	
	 printf("%d\t%f\n",x[i],y[i]);
	
      }; */

//       xvgfile = "rmsf.xvg";
       
//       i = read_xvg_files(xvgfile,&sizeof_x,x,y);
       
//       printf("%s\n","here");
       
       maxi_king = 0;
       
       for(i=1;i<=sizerms;i++)
       {
    	   if(y[i] > maxi_king) maxi_king = y[i];
	   
//	   printf("%f\t%s\n",maxi_king,"y");
       };
       
       class_rank = maxi_king/10; // 10 classes , division of max_rms / 10 ,
       
       for(k=0;k<=10;k++){
	
	   kf = (float)k;
	 
	   test_rank[k] = (float)class_rank*kf; // each class = part of max_rms / 10 * number_of_class
	 
//	   printf("%f\t%s\n",test_rank[k],"testrank");
      };

//      return 0;
//        printf("%s\n","here 2 ");     
      
       for(k=1;k<=10;k++){
      
	 o = 0;
	 
	 empty_class[k] = 1;	 //empty default set
         
	 for(i=1;i <= sizerms; i++) {
	 
	             if(y[i] >= test_rank[k-1] && y[i] <= test_rank[k]){
		
		       o ++;		       

		        class_res[k][o] = x[i];  // test if res i, with resnumber x is in class k. class res = are residues of class with number o
			       
			  rms_res[k][o] = y[i];
			  
			  class_num[k]  = o;
			  
			  empty_class[k] = 0;  //filled
	//		  printf("%d\t%s\n",o,"class num");
		    };
	 
          };
	  
	  
       };	  
        
    /*   for(k=1;k<=10;k++){
	
	 printf("%d\t%s\n",class_num[k],"class num");
	 
	 if(empty_class[k] == 1) printf("%s\n","empty");
	 if(empty_class[k] == 0) printf("%s\n","full");	 
	 
      };*/
       
       // remove empty classes
       
       o = 1;
       
       n = 0;
       
       for(k=1;k<=10;k++){
	 
	   if(empty_class[k] == 0){ // if class is filled
	    
	     for(o=1;o<=class_num[k];o++){
	     
	        if(o == 1){         
		  
		  n ++;             // new counter for class_pupils n
		  
		};   
	       
	        class_pupils[n][o]  = class_res[k][o];    

	        pupil_rms[n][o]     = rms_res[k][o]; 
	        
	        number_of_pupils[n] = o;  // number of pupils in each class n

//	        printf("%d\t%s\n",o,"here do while");		        
	   
	     };
	     
	  };
	 
      };

       
      num_class_real_f = (float)n ; // number of filled classes

     if(num_class_real_f < eventmax) // if rmsd distribution too flat. AT minimum it has to be as large as eventmax
// note eventmax is set to 5 !
        {

            // rigid random selection if num_class_real < eventmax

         for(n=1;n<=10;n++)
             {

             for(o=1;o <= eventmax ; o++) // now each class is filled
                {

                  random_number(&ran);     

    // obsolete and wrong             class_pupils[n][o] =  res_number[(int)roundf(numat_f*ran)];  // random fill
                 class_pupils[n][o] = range_checker_A + (int)roundf(ran*(range_checker_B - range_checker_A)); // random-selection within range A to B.  


                 number_of_pupils[n] =  eventmax ; // setting number of residues in KMC to eventmax // pupil_rms already set

               }; // end of rigid random selection

            }; // end of loop over classes

          

            outkmc = fopen("OUTPUT","a");

                fprintf(outkmc,"%f\t%s\n",num_class_real_f,"lead to rigid distribution to 10 classes");

           fclose(outkmc);

           num_class_real_f = 10; // reset of number of classes to 10


        };    // end of if condition 'eventmax' 

      // ran selection
      
       random_number(&ran);
       
      selected_class = roundf(num_class_real_f* ran) ;   
  
       printf("%d\t%s\n",selected_class,"selected class");
   
       if(selected_class == 0) selected_class = 1;    

       outkmc = fopen("OUTPUT","a"); 
       
                fprintf(outkmc,"%d\t%s\t%f\t%s\n",selected_class,"selected class of",num_class_real_f,"classes");
		
	fclose(outkmc);       

       
       fp = fopen("searchlist.ndx","w");
       
       for(i=1;i<=number_of_pupils[selected_class];i++)
       {
	   fprintf(fp,"%d\t%f\n",class_pupils[selected_class][i],pupil_rms[selected_class][i]); // let selected class go on KMC-trip !
       };
      
       fclose(fp);

       
}




int mdprod(int argc,char *argv[],int mdtime,float Temp,int imp,char *grofile)
{
        int i;
	
	char *topfile;
	 
	char *trajfile,*final_file,*ener_file,*mdpfile;
	
	FILE *fp5;

	
	i = write_mdp(mdtime,Temp,imp);
  
//        system("grompp -f brown.mdp -p start.top -c start.gro -maxwarn 3 -o run.tpr");

//       system("ibrun mdrun_mpi -v -s run.tpr -c start.gro -o md_traj1.trr -e md_ener.edr -pd ");
	  	

//       system("grompp -f langevin.mdp -p start.top -c start.gro -maxwarn 3 -o run.tpr");

//        system("ibrun mdrun_mpi -v -s run.tpr -c start.gro -o md_traj2.trr -e md_ener.edr -pd ");

	
        system("~/GRO453/bin/grompp -f fullmd_sol.mdp -p start.top -c start.gro -maxwarn 3 -o run.tpr "); 	
	

	system("mpirun -np 6 ~/gromacs-4.5.5/src/kernel/mdrun -v -s run.tpr -c minimized3.gro -o md_traj.trr -e md_ener.edr -pd ");
	
/*	fp5 = fopen("mdcat.sh","w");
	
	fprintf(fp5,"%s\n","#!/bin/bash");
	fprintf(fp5,"%s\n","trjcat -f md_traj1.trr md_traj2.trr md_traj3.trr -settime -o md_traj.trr << EOF");
	fprintf(fp5,"%s\n","0");
	fprintf(fp5,"%d\n",mdtime*3);
        fprintf(fp5,"%d\n",mdtime*5);
	fprintf(fp5,"%s\n","EOF");
	
	fclose(fp5);
	
	system("chmod 744 mdcat.sh");
	
	system("./mdcat.sh");
*/	
  
}

int write_mdp(int mdtime, float Temp,int imp)
{
        FILE *fp2;
	
        FILE *fp = fopen("fullmd_sol.mdp","w");
	
	fprintf(fp,"%s\n","title                    = implicit prod.");
	fprintf(fp,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp,"%s\n","include 		 = -I../top");
        fprintf(fp,"%s\n","define                = -DPOSRE");
        fprintf(fp,"%s\n","integrator		 = md");
	fprintf(fp,"%s\n","dt			 = 0.001");
	fprintf(fp,"%s %d\n","nsteps  		 = ",mdtime);
	fprintf(fp,"%s\n","comm_mode         =  linear ");
	fprintf(fp,"%s\n","nstxout 		 = 5000 ");
	fprintf(fp,"%s\n","nstvout 		 = 5000 ");
	fprintf(fp,"%s\n","nstlog  		 = 5000 ");
	fprintf(fp,"%s\n","nstenergy		 = 5000 ");
	fprintf(fp,"%s\n","nstxtcout		 = 5000 ");	
	fprintf(fp,"%s\n","xtc_grps		 = protein");
	fprintf(fp,"%s\n","energygrps		 = protein ");	
	fprintf(fp,"%s\n","nstlist 		 = 2");	
	fprintf(fp,"%s\n","ns_type 		 = grid");		
	fprintf(fp,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp,"%s\n","coulombtype		 = pme ");		
	fprintf(fp,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp,"%s\n","vdwtype           = shift ");	
	fprintf(fp,"%s\n","rvdw			 = 1.0");	
	fprintf(fp,"%s\n","pbc                      = xyz");	
	fprintf(fp,"%s\n","tcoupl  		 = v-rescale");	
	fprintf(fp,"%s\n","tc-grps 		 = system");	
	fprintf(fp,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp,"%s %f \n","ref_t		 = ",Temp);	
	fprintf(fp,"%s\n","Pcoupl  		 = no");
	fprintf(fp,"%s\n","pcoupltype = isotropic");
	fprintf(fp,"%s\n","compressibility = 4.5E-5");
	fprintf(fp,"%s\n","ref_p = 1.0");
	fprintf(fp,"%s\n","tau_p = 1.0");
	fprintf(fp,"%s\n","gen_vel 		 = yes");	
	fprintf(fp,"%s %f \n","gen_temp		 = ",Temp);	
        fprintf(fp,"%s\n","constraints         = none");
	fclose(fp);
	
        fp = fopen("brown.mdp","w");
	
	fprintf(fp,"%s\n","title                    = implicit prod.");
	fprintf(fp,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp,"%s\n","include 		 = -I../top");
        fprintf(fp,"%s\n","define                = -DPOSRE");
        fprintf(fp,"%s\n","integrator		 = bd");
	fprintf(fp,"%s\n","bd-fric              = 5000");
	fprintf(fp,"%s\n","dt			 = 0.003");
	fprintf(fp,"%s %d\n","nsteps  		 = ",mdtime);
	fprintf(fp,"%s\n","comm_mode         =  linear ");
	fprintf(fp,"%s\n","nstxout 		 = 5000 ");
	fprintf(fp,"%s\n","nstvout 		 = 5000 ");
	fprintf(fp,"%s\n","nstlog  		 = 5000 ");
	fprintf(fp,"%s\n","nstenergy		 = 5000 ");
	fprintf(fp,"%s\n","nstxtcout		 = 5000 ");	
	fprintf(fp,"%s\n","xtc_grps		 = protein");
	fprintf(fp,"%s\n","energygrps		 = protein ");	
	fprintf(fp,"%s\n","nstlist 		 = 2");	
	fprintf(fp,"%s\n","ns_type 		 = grid");		
	fprintf(fp,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp,"%s\n","coulombtype		 = pme ");		
	fprintf(fp,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp,"%s\n","vdwtype           = shift ");	
	fprintf(fp,"%s\n","rvdw			 = 1.0");	
	fprintf(fp,"%s\n","pbc                      = xyz");	
	fprintf(fp,"%s\n","tcoupl  		 = v-rescale");	
	fprintf(fp,"%s\n","tc-grps 		 = system");	
	fprintf(fp,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp,"%s %f \n","ref_t		 = ",Temp);	
	fprintf(fp,"%s\n","Pcoupl  		 = no");
	fprintf(fp,"%s\n","pcoupltype = isotropic");
	fprintf(fp,"%s\n","compressibility = 4.5E-5");
	fprintf(fp,"%s\n","ref_p = 1.0");
	fprintf(fp,"%s\n","tau_p = 1.0");
	fprintf(fp,"%s\n","gen_vel 		 = yes");	
	fprintf(fp,"%s %f \n","gen_temp		 = ",Temp);	
        fprintf(fp,"%s\n","constraints         = H-bonds");
	fclose(fp);	
	
        fp = fopen("langevin.mdp","w");
	
	fprintf(fp,"%s\n","title                    = implicit prod.");
	fprintf(fp,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp,"%s\n","include 		 = -I../top");
        fprintf(fp,"%s\n","define                = -DPOSRE");
        fprintf(fp,"%s\n","integrator		 = bd");
	fprintf(fp,"%s\n","bd-fric              = 1000");
	fprintf(fp,"%s\n","dt			 = 0.001");
	fprintf(fp,"%s %d\n","nsteps  		 = ",mdtime);
	fprintf(fp,"%s\n","comm_mode         =  linear ");
	fprintf(fp,"%s\n","nstxout 		 = 5000 ");
	fprintf(fp,"%s\n","nstvout 		 = 5000 ");
	fprintf(fp,"%s\n","nstlog  		 = 5000 ");
	fprintf(fp,"%s\n","nstenergy		 = 5000 ");
	fprintf(fp,"%s\n","nstxtcout		 = 5000 ");	
	fprintf(fp,"%s\n","xtc_grps		 = protein");
	fprintf(fp,"%s\n","energygrps		 = protein ");	
	fprintf(fp,"%s\n","nstlist 		 = 2");	
	fprintf(fp,"%s\n","ns_type 		 = grid");		
	fprintf(fp,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp,"%s\n","coulombtype		 = pme ");		
	fprintf(fp,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp,"%s\n","vdwtype           = shift ");	
	fprintf(fp,"%s\n","rvdw			 = 1.0 ");	
	fprintf(fp,"%s\n","pbc                      = xyz");	
	fprintf(fp,"%s\n","tcoupl  		 = v-rescale");	
	fprintf(fp,"%s\n","tc-grps 		 = system");	
	fprintf(fp,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp,"%s %f \n","ref_t		 = ",Temp);	
	fprintf(fp,"%s\n","Pcoupl  		 = no");
	fprintf(fp,"%s\n","pcoupltype = isotropic");
	fprintf(fp,"%s\n","compressibility = 4.5E-5");
	fprintf(fp,"%s\n","ref_p = 1.0");
	fprintf(fp,"%s\n","tau_p = 1.0");
	fprintf(fp,"%s\n","gen_vel 		 = yes");	
	fprintf(fp,"%s %f \n","gen_temp		 = ",Temp);	
	fclose(fp);	

	fp2 = fopen("minim.mdp","w");
	
	fprintf(fp2,"%s\n","title                    = implicit minim.");
	fprintf(fp2,"%s\n","cpp			 = /lib/cpp");
	fprintf(fp2,"%s\n","include 		 = -I../top");
        fprintf(fp2,"%s\n","define               = -DPOSRE");
        fprintf(fp2,"%s\n","integrator		 = l-bfgs");
	fprintf(fp2,"%s\n","emtol                = 2000");
	fprintf(fp2,"%s\n","emstep               = 0.02");
	fprintf(fp2,"%s\n","nsteps               = 200");
	fprintf(fp2,"%s\n","dt			 = 0.001");
	fprintf(fp2,"%s\n","comm_mode            = linear ");
	fprintf(fp2,"%s\n","nstxout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstvout 		 = 1500 ");
	fprintf(fp2,"%s\n","nstlog  		 = 1500 ");
	fprintf(fp2,"%s\n","nstenergy		 = 1500 ");
	fprintf(fp2,"%s\n","nstxtcout		 = 1500 ");	
	fprintf(fp2,"%s\n","xtc_grps		 = protein");
	fprintf(fp2,"%s\n","energygrps		 = protein ");	
	fprintf(fp2,"%s\n","nstlist 		 = 10");	
	fprintf(fp2,"%s\n","ns_type 		 = grid");		
	fprintf(fp2,"%s\n","rlist		 = 1.0 ");	
	fprintf(fp2,"%s\n","coulombtype		 = pme");		
	fprintf(fp2,"%s\n","rcoulomb		 = 1.0");
	fprintf(fp2,"%s\n","vdwtype           = shift ");	
	fprintf(fp2,"%s\n","rvdw			 = 1.0");	
	fprintf(fp2,"%s\n","pbc                  = xyz");	
	fprintf(fp2,"%s\n","tcoupl  		 = no");	
	fprintf(fp2,"%s\n","tc-grps 		 = system");	
	fprintf(fp2,"%s\n","tau_t	        =  1.0 ");	
	fprintf(fp2,"%s %f \n","ref_t		 = ",Temp);	
	fprintf(fp2,"%s\n","Pcoupl  		 = no");	
	fprintf(fp2,"%s\n","gen_vel 		 = no");	
	fprintf(fp2,"%s %f \n","gen_temp		 = ",Temp);	
        fprintf(fp2,"%s\n","constraints         = none");
	fclose(fp2);
	
return 1;	
}

int replex(int argc, char *argv[],float T1,float T2,float T3,float T4,float T5,float *sim_temp,int mdtime,float Temp_orig,int imp)
{

  int i,k,u;

  float R = 8.314E-3;
  
  int Temp_number;
  
  float Temp_array[10];
  
  float Temp_simulated[2];
  
  float  rep_energy[10];
  
  float  en_aver;
  
  float   x;
  
  int seed,y,select_term;
  
  int mdtime2,mpierr;
  
  char *trajfile,*final_file,*ener_file,*grofile,*topfile;
  
  char *mdpfile;
  
  mdtime2   = 10000;

  argc = 1;

 char *xvgfile;
float aver,sizeof_x2;
int sizeof_x;
int x2[100000];
float y2[100000];	
 
  FILE *fp;
 
         i = write_mdp(mdtime2,Temp_orig,imp);

  
    system("~/GRO453/bin/grompp -f minim.mdp -p start.top -c start.gro  -maxwarn 3 -o run.tpr"); 	
	

 	system("~/GRO453/bin/mdrun -v -s run.tpr -c start.gro -o trj.trr -e ener.edr ");
  
  grofile   = "start";
  trajfile  = "trj";
  ener_file  = "md_ener";
  topfile   = "start";
  final_file = "replex";
  mdpfile   = "fullmd_sol";
  select_term = 11;
  
  for(i=1;i <= 5;i++) {
    
    if(i == 1) 
      
      Temp_array[i] = T1;
    
    if(i == 2) 
      
      Temp_array[i] = T2;    

    if(i == 3) 
      
      Temp_array[i] = T3;

    if(i == 4) 
      
      Temp_array[i] = T4;    
    
    if(i == 5) 
      
      Temp_array[i] = T5;
    
  };
    

 for(i=1;i <= 5;i++) {

   if(Temp_array[i] == Temp_orig) 

      Temp_number = i;
 }; 


   i = random_number(&x);
   

     if(x <= 0.5){
       
       for(y=1;y <= 2; y++) 
	 
       {

	 
	 i = write_mdp(mdtime2,Temp_array[Temp_number+y-1],imp);
	 
	
        system("~/GRO453/bin/grompp -f fullmd_sol.mdp -p start.top -c start.gro  -maxwarn 3 -o run.tpr"); 	
	
	system("~/GRO453/bin/mdrun -v -s run.tpr -c replex.gro -o trj.trr -e ener.edr -pd");	 
	 
	 select_term = 11;
   
/*
	 argc = 5;
	 argv[0] = "g_energy";
	 argv[1] = "-f";
	 argv[2] = "ener.edr";
	 argv[3] = "-o";
	 argv[4] = "1.xvg";
	 
	 i = gmx_energy(argc,argv,ener_file,select_term,&en_aver);*/

    fp = fopen("en.sh","w");
    fprintf(fp,"%s\n","#!/bin/bash");
    fprintf(fp,"%s\n","~/GRO453/bin/g_energy -f md_ener.edr -o t.xvg << EOF");
    fprintf(fp,"%s\n","13");
    fprintf(fp,"%s\n","0");
    fprintf(fp,"%s\n","EOF");
    fclose(fp);

    system("chmod 744 ./en.sh");

    system("./en.sh");


    xvgfile = "t.xvg";    
 
    read_xvg_files(xvgfile,&sizeof_x,x2,y2);

    aver = 0;

    for(k=1;k<=sizeof_x;k++){

                 aver = aver + y2[k];

    };

    sizeof_x2 = (float)sizeof_x;

    en_aver = aver / sizeof_x2;
		 
	 rep_energy[y] = en_aver;
	 
	 Temp_simulated[y] = Temp_array[Temp_number+y-1];
	 
	 

	 };
	 
       };
 
     
     if(x > 0.5){
       
       
       for(y=1;y <= 2; y++) 
	 
       {

	
	 i = write_mdp(mdtime2,Temp_array[Temp_number-y+1],imp);
	 
    system("~/GRO453/bin/grompp -f fullmd_sol.mdp -p start.top -c start.gro  -maxwarn 3 -o run.tpr"); 	

	system("~/GRO453/bin/mdrun -v -s run.tpr -c replex.gro -o trj.trr -e ener.edr -pd");	

	 /*
	 argc = 5;
	 argv[0] = "g_energy";
	 argv[1] = "-f";
	 argv[2] = "ener.edr";
	 argv[3] = "-o";
	 argv[4] = "1.xvg";
	 
	 i = gmx_energy(argc,argv,ener_file,select_term,&en_aver);*/

    fp = fopen("en.sh","w");
    fprintf(fp,"%s\n","#!/bin/bash");
    fprintf(fp,"%s\n","~/GRO453/bin/g_energy -f md_ener.edr -o t.xvg << EOF");
    fprintf(fp,"%s\n","13");
    fprintf(fp,"%s\n","0");
    fprintf(fp,"%s\n","EOF");
    fclose(fp);

    system("chmod 744 ./en.sh");

    system("./en.sh");


    xvgfile = "t.xvg";    
 
    read_xvg_files(xvgfile,&sizeof_x,x2,y2);

    aver = 0;

    for(k=1;k<=sizeof_x;k++){

                 aver = aver + y2[k];

    };

    sizeof_x2 = (float)sizeof_x;

    en_aver = aver / sizeof_x2;
	 
	 rep_energy[y] = en_aver;
	 
	 Temp_simulated[y] = Temp_array[Temp_number-y+1];
	 
	 
	 
       };
	 
      };
      
        
       
      if(Temp_number == 1){

	for(y=1;y <= 2; y++) 
	 
       {

	
	 i = write_mdp(mdtime2,Temp_array[Temp_number+y-1],imp);
	 
     system("~/GRO453/bin/grompp -f fullmd_sol.mdp -p start.top -c start.gro  -maxwarn 3 -o run.tpr"); 
       
	system("~/GRO453/bin/mdrun -v -s run.tpr -c replex.gro -o trj.trr -e ener.edr -pd");	 
	
	 select_term = 11;

	   
	/* argc = 5;
	 argv[0] = "g_energy";
	 argv[1] = "-f";
	 argv[2] = "ener.edr";
	 argv[3] = "-o";
	 argv[4] = "1.xvg";	
	 
	 i = gmx_energy(argc,argv,ener_file,select_term,&en_aver);*/
    fp = fopen("en.sh","w");
    fprintf(fp,"%s\n","#!/bin/bash");
    fprintf(fp,"%s\n","~/GRO453/bin/g_energy -f md_ener.edr -o t.xvg << EOF");
    fprintf(fp,"%s\n","13");
    fprintf(fp,"%s\n","0");
    fprintf(fp,"%s\n","EOF");
    fclose(fp);

    system("chmod 744 ./en.sh");

    system("./en.sh");

    xvgfile = "t.xvg";    
 
    read_xvg_files(xvgfile,&sizeof_x,x2,y2);

    aver = 0;

    for(k=1;k<=sizeof_x;k++){

                 aver = aver + y2[k];

    };

    sizeof_x2 = (float)sizeof_x;

    en_aver = aver / sizeof_x2;
	 
	 rep_energy[y] = en_aver;
	 
	 Temp_simulated[y] = Temp_array[Temp_number+y-1];
	 
         };
	 
       };


      if(Temp_number == 5){

	for(y=1;y <= 2; y++) 
	 
       {

	 
	 i = write_mdp(mdtime2,Temp_array[Temp_number-y+1],imp);
	 

    system("~/GRO453/bin/grompp -f fullmd_sol.mdp -p start.top -c start.gro  -maxwarn 3 -o run.tpr"); 	
	
	system("~/GRO453/bin/mdrun -v -s run.tpr -c replex.gro -o trj.trr -e ener.edr -pd");	
	 
	 select_term = 11;

	/*   
	 argc = 5;
	 argv[0] = "g_energy";
	 argv[1] = "-f";
	 argv[2] = "ener.edr";
	 argv[3] = "-o";
	 argv[4] = "1.xvg";	   
	 
	 i = gmx_energy(argc,argv,ener_file,select_term,&en_aver);*/

    fp = fopen("en.sh","w");
    fprintf(fp,"%s\n","#!/bin/bash");
    fprintf(fp,"%s\n","~/GRO453/bin/g_energy -f md_ener.edr -o t.xvg << EOF");
    fprintf(fp,"%s\n","13");
    fprintf(fp,"%s\n","0");
    fprintf(fp,"%s\n","EOF");
    fclose(fp);

    system("chmod 744 ./en.sh");

    system("./en.sh");


    xvgfile = "t.xvg";    
 
    read_xvg_files(xvgfile,&sizeof_x,x2,y2);

    aver = 0;

    for(k=1;k<=sizeof_x;k++){

                 aver = aver + y2[k];

    };

    sizeof_x2 = (float)sizeof_x;

    en_aver = aver / sizeof_x2;	 
	 rep_energy[y] = en_aver;
	 
	 Temp_simulated[y] = Temp_array[Temp_number-y+1];
	 
	 };
	 
       };
       
      
       
       i = 1;

       if( rep_energy[2] <= rep_energy[1]) 

                   *sim_temp = Temp_simulated[2]
        ;

        if( rep_energy[2] >= rep_energy[1])

                   i = random_number(&x); 

                   if(x >= exp(((1/(R*Temp_simulated[1]))-(1/(R*Temp_simulated[2])))*(rep_energy[1]-rep_energy[2]))) 

                      *sim_temp = Temp_simulated[2];

                   else

                      *sim_temp = Temp_simulated[1]

                   ;

         ;
	 
system("rm -f ./#*");



}

int random_number(float *x)

{
   
/*   int seed;
  
   seed = time(NULL);
   srand(seed);  
   
      
   srand( (unsigned)time( NULL ) ); */ 
   *x = (float) rand()  / (float) RAND_MAX;
   
   if(*x <= 0.09) *x = *x * 10;
   
   if(*x <= 0.009) *x = *x * 100;
   
   if(*x <= 0.0009) *x = *x * 1000;
   
   if(*x <= 0.00009) *x = *x * 10000;
   
   if(*x <= 0.000009) *x = *x * 100000;
   
   if(*x <= 0.0000009) *x = *x * 1000000;   
   
      
}
