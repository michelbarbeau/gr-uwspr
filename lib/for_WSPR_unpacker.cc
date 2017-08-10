   // In ".h"
   FILE *fp_fftwf_wisdom_file, *fall_wspr, *fwsprd, *fhash, *ftimer;
   strcpy(all_fname,".");
   char all_fname[200],timer_fname[200];

 
   // In ".cc"


   // INITIALIZATION

   strcpy(all_fname,".");
   strcpy(spots_fname,".");
   if(data_dir != NULL) {
      strcpy(all_fname,data_dir);
      strcpy(spots_fname,data_dir);
      strcpy(timer_fname,data_dir);
   }
   strncat(all_fname,"/ALL_WSPR.TXT",20);
   fall_wspr=fopen(all_fname,"a");
   strncat(spots_fname,"/wspr_spots.txt",20);
   fwsprd=fopen(spots_fname,"w");

   strncat(timer_fname,"/wspr_timer.out",20);
   if((ftimer=fopen(timer_fname,"r"))) {
      //Accumulate timing data
      nr=fscanf(ftimer,"%f %f %f %f %f %f %f",
                &treadwav,&tcandidates,&tsync0,&tsync1,&tsync2,&tfano,&ttotal);
      fclose(ftimer);
   }
   ftimer=fopen(timer_fname,"w");

   // DESCTRUCTOR
   fclose(fall_wspr);
   fclose(fwsprd);
   fclose(ftimer);  


   // AFTER MAIN LOOP


   // sort the result in order of increasing frequency
   struct result temp;
   for (j = 1; j <= uniques - 1; j++) {
      for (k = 0; k < uniques - j; k++) {
         if (decodes[k].freq > decodes[k+1].freq) {
            temp = decodes[k];
            decodes[k]=decodes[k+1];;
            decodes[k+1] = temp;
         }
      }
   }

   // print decoded data, in files
   // stdout
   // fall_wspr = "/ALL_WSPR.TXT"
   // fwsprd = "/wspr_spots.txt"
   for (i=0; i<uniques; i++) {
      printf("%4s %3.0f %4.1f %10.6f %2d  %-s \n",
             decodes[i].time, decodes[i].snr,decodes[i].dt, decodes[i].freq,
             (int)decodes[i].drift, decodes[i].message);
      fprintf(fall_wspr,
              "%6s %4s %3d %3.0f %5.2f %11.7f  %-22s %2d %5u %4d\n",
              decodes[i].date, decodes[i].time, (int)(10*decodes[i].sync),
              decodes[i].snr, decodes[i].dt, decodes[i].freq,
              decodes[i].message, (int)decodes[i].drift, decodes[i].cycles/81,
              decodes[i].jitter);
      fprintf(fwsprd,
              "%6s %4s %3d %3.0f %4.1f %10.6f  %-22s %2d %5u %4d\n",
              decodes[i].date, decodes[i].time, (int)(10*decodes[i].sync),
              decodes[i].snr, decodes[i].dt, decodes[i].freq,
              decodes[i].message, (int)decodes[i].drift, decodes[i].cycles/81,
              decodes[i].jitter);

   }


   // print time stats (in file ftimer="wspr_timer.out")
   ttotal += (float)(clock()-t00)/CLOCKS_PER_SEC;

   fprintf(ftimer,"%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n\n",
           treadwav,tcandidates,tsync0,tsync1,tsync2,tfano,ttotal);

   fprintf(ftimer,"Code segment        Seconds   Frac\n");
   fprintf(ftimer,"-----------------------------------\n");
   fprintf(ftimer,"readwavfile        %7.2f %7.2f\n",treadwav,treadwav/ttotal);
   fprintf(ftimer,"Coarse DT f0 f1    %7.2f %7.2f\n",tcandidates,
           tcandidates/ttotal);
   fprintf(ftimer,"sync_and_demod(0)  %7.2f %7.2f\n",tsync0,tsync0/ttotal);
   fprintf(ftimer,"sync_and_demod(1)  %7.2f %7.2f\n",tsync1,tsync1/ttotal);
   fprintf(ftimer,"sync_and_demod(2)  %7.2f %7.2f\n",tsync2,tsync2/ttotal);
   fprintf(ftimer,"Stack/Fano decoder %7.2f %7.2f\n",tfano,tfano/ttotal);
   fprintf(ftimer,"-----------------------------------\n");
   fprintf(ftimer,"Total              %7.2f %7.2f\n",ttotal,1.0);
