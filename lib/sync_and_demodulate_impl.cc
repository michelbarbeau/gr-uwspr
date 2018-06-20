/* -*- c++ -*- */
/*
 * Copyright 2018 Michel Barbeau, Carleton University.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "sync_and_demodulate_impl.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <stdbool.h>
#include <volk/volk.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
#include <assert.h>

#ifndef DEBUG
#define DEBUG 1 // set debug mode
#endif
#include "debugmacro.h"

using namespace std;

namespace gr {
 namespace uwspr {

   sync_and_demodulate::sptr
   sync_and_demodulate::make(int fs, int fl, int spb, int maxdrift,
     int maxfreqs, int cf)
   {
     return gnuradio::get_initial_sptr
      (new sync_and_demodulate_impl(fs, fl, spb, maxdrift, maxfreqs, cf));
   }

/*
 * The private constructor
 */
   sync_and_demodulate_impl::sync_and_demodulate_impl(
     int fs, int fl, int spb, int maxdrift, int maxfreqs, int cf)
      : gr::block("sync_and_demodulate",
        gr::io_signature::make(0, 0, 0),
        gr::io_signature::make(0, 0, 0))
   {
     // create and register the asynchronous input port
     in_port = pmt::mp("in");
     message_port_register_in(in_port);
     set_msg_handler(in_port, boost::bind(
       &sync_and_demodulate_impl::demodulate, this,_1));
     // create and register the asynchronous output port
     out_port = pmt::mp("out");
     message_port_register_out(out_port);
     // init sampling frequency
     this->fs = fs;
     // init frame length, in samples
     this->fl = fl;
     // init number of samples by baud (symbol)
     this->spb = spb;
     // init number of samples by baud (symbol)
     this->maxdrift = maxdrift;
     // maximum number of frequencies consider during the search
     this->maxfreqs = maxfreqs;
     // carrier frequency
     this->carrierfrequency = cf;
     // init array of candidate frequencies
     candidates = new candidate_t[maxfreqs];
     // init delta frequency
     df=375.0/256.0/2; // or df=775/256
     npoints = 45000;
     nffts=4*floor(npoints/512)-1;
     // metric table initialization
#include "./metric_tables.c"
     bias=0.45;    //Fano metric bias (used for both Fano and stack algorithms)
     for(int i=0; i<256; i++) {
       mettab[0][i]=round( 10*(metric_tables[2][i]-bias) );
       mettab[1][i]=round( 10*(metric_tables[2][255-i]-bias) );
     }
     tfano=0.0; treadwav=0.0; tcandidates=0.0; tsync0=0.0;
     tsync1=0.0; tsync2=0.0; ttotal=0.0;
     nbits=81;
     symbols=(unsigned char *)malloc(sizeof(char)*nbits*2);
     decdata=(unsigned char *)malloc(sizeof(char)*11);
     channel_symbols=(unsigned char *)malloc(sizeof(char)*nbits*2);
     // open message log file
     msglogfile = fopen("messagelog.txt","a");
     if (msglogfile)
     {
       fprintf(stderr, "Messages logged in file messagelog.txt\n");
     } else {
     fprintf(stderr, "Error opening message log file!\n");
     }
     time(&start);
     info = localtime(&start);
     fprintf(msglogfile, "Start time: %s\n", asctime(info));
     fflush (msglogfile);
     // init frame count
     framecount = 0;
   }

   /*
    * Virtual destructor.
    */
   sync_and_demodulate_impl::~sync_and_demodulate_impl()
   {
     fclose(msglogfile);
   }

/* WSPR uses the Layland-Lushbaugh code
 * Nonsystematic, non-quick look-in, dmin=?, dfree=?
 */
#define POLY1 0xf2d05351
#define POLY2 0xe4613c47

/* Convolutional encoder macro. Takes the encoder state, generates
 * a rate 1/2 symbol pair and stores it in 'sym'. The symbol generated from
 * POLY1 goes into the 2-bit of sym, and the symbol generated from POLY2
 * goes into the 1-bit.
 * Copyright 1994, Phil Karn, KA9Q
 */
#define ENCODE(sym,encstate){ \
 unsigned long _tmp; \
 \
 _tmp = (encstate) & POLY1; \
 _tmp ^= _tmp >> 16; \
 (sym) = Partab[(_tmp ^ (_tmp >> 8)) & 0xff] << 1; \
 _tmp = (encstate) & POLY2; \
 _tmp ^= _tmp >> 16; \
 (sym) |= Partab[(_tmp ^ (_tmp >> 8)) & 0xff]; \
}

   /* Convolutionally encode a packet. The input data bytes are read
    * high bit first and the encoded packet is written into 'symbols',
    * one symbol per byte. The first symbol is generated from POLY1,
    * the second from POLY2.
    * Storing only one symbol per byte uses more space, but it is faster
    * and easier than trying to pack them more compactly.
    */
   int sync_and_demodulate_impl::encode(
        unsigned char *symbols,     // Output buffer, 2*nbytes*8
        unsigned char *data,     // Input buffer, nbytes
        unsigned int nbytes)     // Number of bytes in data
   {
     unsigned long encstate;
     int sym;
     int i;
     encstate = 0;
     while(nbytes-- != 0) {
       for(i=7; i>=0; i--) {
         encstate = (encstate << 1) | ((*data >> i) & 1);
         ENCODE(sym,encstate);
         *symbols++ = sym >> 1;
         *symbols++ = sym & 1;
       }
       data++;
     }
     return 0;
   }

   /*
    * Soft decision Fano sequential decoder for K=32 r=1/2 convolutional code
    * Copyright 1994, Phil Karn, KA9Q
    * Decode packet with the Fano algorithm.
    * Return 0 on success, -1 on timeout
    * Decode packet with the Fano algorithm.
    * Return 0 on success, -1 on timeout
    */
   int sync_and_demodulate_impl::fano(
     unsigned int  *metric, // Final path metric (returned value)
     unsigned int  *cycles, // Cycle count (returned value)
     unsigned int  *maxnp, // Progress before timeout (returned value)
     unsigned char *data, // Decoded output data
     unsigned char *symbols, // Raw deinterleaved input symbols
     unsigned int nbits, // Number of output bits
     int mettab[2][256], // Metric table, [sent sym][rx symbol]
     int delta, // Threshold adjust parameter
     unsigned int maxcycles) // Decoding timeout in cycles per bit
   {
     struct node *nodes; // First node
     struct node *np; // Current node
     struct node *lastnode; // Last node
     struct node *tail; // First node of tail
     int t;       // Threshold
     int m0,m1;
     int ngamma;
     unsigned int lsym;
     unsigned int i;
     if((nodes = (struct node *)malloc((nbits+1)*sizeof(struct node))) == NULL)
     {
       printf("malloc failed\n");
       return 0;
     }
     lastnode = &nodes[nbits-1];
     tail = &nodes[nbits-31];
     *maxnp = 0;
     /* Compute all possible branch metrics for each symbol pair
      * This is the only place we actually look at the raw input symbols
      */
     for(np=nodes; np <= lastnode; np++)
     {
       np->metrics[0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
       np->metrics[1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
       np->metrics[2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
       np->metrics[3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
       symbols += 2;
     }
     np = nodes;
     np->encstate = 0;
     // Compute and sort branch metrics from root node */
     ENCODE(lsym,np->encstate);     // 0-branch (LSB is 0)
     m0 = np->metrics[lsym];
     /* Now do the 1-branch. To save another ENCODE call here and
      * inside the loop, we assume that both polynomials are odd,
      * providing complementary pairs of branch symbols.
      * This code should be modified if a systematic code were used.
      */
     m1 = np->metrics[3^lsym];
     if(m0 > m1) {
       np->tm[0] = m0; // 0-branch has better metric
       np->tm[1] = m1;
     } else {
       np->tm[0] = m1; // 1-branch is better
       np->tm[1] = m0;
       np->encstate++; // Set low bit
     }
     np->i = 0; // Start with best branch
     maxcycles *= nbits;
     np->gamma = t = 0;
     // Start the Fano decoder
     for(i=1; i <= maxcycles; i++) {
       if((int)(np-nodes) > (int)*maxnp) *maxnp=(int)(np-nodes);
#ifdef  debug
       printf("k=%ld, g=%ld, t=%d, m[%d]=%d, maxnp=%d, encstate=%lx\n",
         np-nodes,np->gamma,t,np->i,np->tm[np->i],*maxnp,np->encstate);
#endif
       // Look forward
       ngamma = np->gamma + np->tm[np->i];
       if(ngamma >= t) {
         if(np->gamma < t + delta) { // Node is acceptable
           /* First time we've visited this node;
            * Tighten threshold.
            *   t += delta * ((ngamma - t)/delta);
            * but the multiply and divide are slower.
            */
           while(ngamma >= t + delta) t += delta;
         }
         np[1].gamma = ngamma; // Move forward
         np[1].encstate = np->encstate << 1;
         if( ++np == (lastnode+1) ) {
           break; // Done!
         }
         /* Compute and sort metrics, starting with the
          * zero branch
          */
         ENCODE(lsym,np->encstate);
         if(np >= tail) {
           /* The tail must be all zeroes, so don't
            * bother computing the 1-branches here.
            */
           np->tm[0] = np->metrics[lsym];
         } else {
           m0 = np->metrics[lsym];
           m1 = np->metrics[3^lsym];
           if(m0 > m1) {
             np->tm[0] = m0; // 0-branch is better
             np->tm[1] = m1;
           } else {
             np->tm[0] = m1; // 1-branch is better
             np->tm[1] = m0;
             np->encstate++; // Set low bit
           }
         }
         np->i = 0; // Start with best branch
         continue;
       }
       // Threshold violated, can't go forward
       for(;; ) { // Look backward
         if(np == nodes || np[-1].gamma < t) {
           /* Can't back up either.
            * Relax threshold and and look
            * forward again to better branch.
            */
           t -= delta;
           if(np->i != 0) {
             np->i = 0;
             np->encstate ^= 1;
           }
           break;
         }
         // Back up
         if(--np < tail && np->i != 1) {
           np->i++; // Search next best branch
           np->encstate ^= 1;
           break;
         } // else keep looking back
       }
     }
     *metric =  np->gamma; // Return the final path metric
     // Copy decoded data to user's buffer
     nbits >>= 3;
     np = &nodes[7];
     while(nbits-- != 0) {
       *data++ = np->encstate;
       np += 8;
     }
     *cycles = i+1;
     free(nodes);
     if(i >= maxcycles) return -1; // Decoder timed out
     return 0; // Successful completion
   }

#include "pr3.h"

   void sync_and_demodulate_impl::sync_and_demodulate(
     candidate_t candidate,
     float *id, float *qd, long np,
     unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
     int *shift1, int lagmin, int lagmax, int lagstep,
     float *drift1, int symfac, float *sync, int mode)
   {
     /***********************************************************************
      * mode = 0: no frequency or drift search. find best time lag.          *
      *        1: no time lag or drift search. find best frequency.          *
      *        2: no frequency or time lag search. calculate soft-decision   *
      *           symbols using passed frequency and shift.                  *
      ************************************************************************/
     // Input signal:
     //    id = in-phase (float array of maxpts=65536)
     //    qd = quarature (float array of maxpts=65536)
     //    np = number of points
     // symbols = (char array of size nbits (81) times 2)
     static float fplast=-10000.0;
     // delta time (dt) and delta (frequency)
     static float dt=1.0/375.0, df=375.0/256.0;
     // baseband frequencies of symbols
     float delta[] = { -df*1.5, -df*0.5, df*0.5, df*1.5 };
     int i, j, k, lag, n;
     float inp[4][162], quad[4][162];
     float p[4]; // power correlated with symbols
     float cmet,totp,syncmax,fac;
     float c[4][256], s[4][256];
     float cdphi, sdphi;
     float ss; // correlation with synchronization soft symbols
     float f0=0.0, fp, fbest=0.0, fsum=0.0, f2sum=0.0, fsymb[162];
     int best_shift = 0, ifreq;
     float t;
     syncmax=-1e30;
     if( mode == 0 ) {ifmin=0; ifmax=0; fstep=0.0; f0=*f1; }
     if( mode == 1 ) {lagmin=*shift1; lagmax=*shift1; f0=*f1; }
     if( mode == 2 ) {lagmin=*shift1; lagmax=*shift1; ifmin=0; ifmax=0; f0=*f1; }
     for(ifreq=ifmin; ifreq<=ifmax; ifreq++) { // frequency search
       f0=*f1+ifreq*fstep;
       for(lag=lagmin; lag<=lagmax; lag=lag+lagstep) { // time search
         ss=0.0; totp=0.0;
         for (i=0; i<162; i++) { // i = frame symbol index
           assert(candidate.m_type==linear||
             candidate.m_type==nonlinear);
           switch(candidate.m_type) {
             // linear drift search
             case linear : {
                fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;
                break;
              }
              // nonlinear drift search
              t = i * 111 / 162; // map symbol index to delay in seconds
              case nonlinear : {
                fp = f0 + slmFrequencyDrift(candidate.m_nonlinear,
                  carrierfrequency, t);
                break;
               }
           }
           // do for each symbol, only if frequency is changed
           if( i==0 || (fp != fplast) ) {
             for (j=0; j<4; j++) {
               // calculate waveform for each sample of symbol
               cdphi=cos(2*M_PI*dt*(fp+delta[j])); // cosine increment of a sample
               sdphi=sin(2*M_PI*dt*(fp+delta[j])); // sine increment of a sample
               c[j][0]=1; s[j][0]=0; // start at phase 0 radian
               for (k=1; k<256; k++) { // k = sample in symbol index
                 // in-phase (cosine)
                 c[j][k]=c[j][k-1]*cdphi - s[j][k-1]*sdphi;
                 // quadrature (sine)
                 s[j][k]=c[j][k-1]*sdphi + s[j][k-1]*cdphi;
               }
               fplast = fp;
             }
           }
           for (j=0; j<4; j++) { // j = defined symbol index
             inp[j][i]=0.0; quad[j][i]=0.0;
             // correlation
             for (k=0; k<256; k++) { // k = sample in symbol index
               n=lag+i*256+k; // n = sample index in time domain representation
               if( (n>0) && (n<np) ) {
                 inp[j][i]=inp[j][i] + id[n]*c[j][k] + qd[n]*s[j][k];
                 quad[j][i]=quad[j][i] - id[n]*s[j][k] + qd[n]*c[j][k];
               }
             }
             // magnitude
             p[j] = sqrt(inp[j][i]*inp[j][i] + quad[j][i]*quad[j][i]);
           }
           totp=totp+p[0]+p[1]+p[2]+p[3]; // total power
           cmet=(p[1]+p[3])-(p[0]+p[2]); // metric
           ss = (pr3[i] == 1) ? ss+cmet : ss-cmet;
           if( mode == 2) { //Compute soft symbols (ss)
             if(pr3[i]==1) {
               // synchro bit is 1
               fsymb[i]=p[3]-p[1]; // difference of powers at data bits 1 and 0
             } else {
               // synchro bit is 0
               fsymb[i]=p[2]-p[0]; // difference of powers at data bits 1 and 0
             }
           }
         }
         ss=ss/totp;
         if( ss > syncmax ) { // Save best parameters
           syncmax=ss;
           best_shift=lag;
           fbest=f0;
         }
       } // lag loop
     } // freq loop
     if( mode <=1 ) { // send best params back to caller
       *sync=syncmax;
       *shift1=best_shift;
       *f1=fbest;
       return;
     }
     if( mode == 2 ) {
       *sync=syncmax;
       for (i=0; i<162; i++) { // normalize the soft symbols
         fsum=fsum+fsymb[i]/162.0;
         f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
       }
       fac=sqrt(f2sum-fsum*fsum);
       for (i=0; i<162; i++) {
         fsymb[i]=symfac*fsymb[i]/fac;
         if( fsymb[i] > 127) fsymb[i]=127.0;
         if( fsymb[i] < -128 ) fsymb[i]=-128.0;
         symbols[i]=fsymb[i] + 128;
       }
       return;
     }
     return;
   }

   int sync_and_demodulate_impl::floatcomp(const void* elem1, const void* elem2)
   {
     if(*(const float*)elem1 < *(const float*)elem2)
        return -1;
     return *(const float*)elem1 > *(const float*)elem2;
   }

   void sync_and_demodulate_impl::deinterleave(unsigned char *sym)
   {
     unsigned char tmp[162];
     unsigned char p, i, j;
     p=0;
     i=0;
     while (p<162) {
       j=((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
       if (j < 162 ) {
         tmp[p]=sym[j];
         p=p+1;
       }
       i=i+1;
     }
     for (i=0; i<162; i++) {
       sym[i]=tmp[i];
     }
   }

   // log a message
   void sync_and_demodulate_impl::logmessage(signed char message[])
   {
     int n;
     const char hex[]={'0','1','2','3','4','5','6','7','8','9','a','b',
       'c','d','e','f'};
       for (int i=0; i<7; i++ ) {
         int n=message[i];
          if(n<0) {
             n = 256+n;
          }
          fprintf(msglogfile,"%c%c", hex[n/16], hex[n%16]);
       }
   }

   // print current time and ellapsed timed
   void sync_and_demodulate_impl::printtime()
   {
     time_t now;
     struct tm *info;
     // handoff work to next block
     time(&now);
     info = localtime(&now);
     fprintf(msglogfile, "Handoff time : %s", asctime(info));
     time_t dtime = difftime(now,start);
     int h = (dtime / 3600) % 24;
     int m = (dtime / 60) % 60;
     int s = dtime % 60;
     fprintf(msglogfile, "Elapsed time: %02d:%02d:%02d\n", h, m, s);
   }

   void sync_and_demodulate_impl::demodulate(pmt::pmt_t msg)
   {
     /* Local variables
      */
     int i,j,k;
     signed char message[]={-9,13,-35,123,57,-39,64,0,0,0,0};
     int delta, maxpts=65536;
     int shift1, lagmin, lagmax, lagstep, ifmin, ifmax, worth_a_try, not_decoded;
     unsigned int metric, cycles, maxnp;
     float df=375.0/256.0/2; // or df=775/256
     float dt=1.0/375.0;
     float fmin=-110, fmax=110;
     float f1, fstep, sync1, drift1;
     // Parameters used for performance-tuning:
     unsigned int maxcycles=10000;             //Decoder timeout limit
     float minsync1=0.10;                      //First sync limit
     float minsync2=0.12;                      //Second sync limit
     int iifac=8;                              //Step size in final DT peakup
     int symfac=50;                            //Soft-symbol normalizing factor
     float minrms=52.0 * (symfac/64.0);       //Final test for plausible decoding
     delta=60;                                //Fano threshold step
     // Get info on candidate frequencies in polymorphic data PDU
     pmt::pmt_t meta(pmt::car(msg)); // CAR is NULL (dictionary)
     pmt::pmt_t tuple(pmt::cdr(msg)); // CDR is a tuple
     // 1st tuple component is signal time domain representation
     pmt::pmt_t vector(pmt::tuple_ref(tuple,0));
     std::complex<double> c; float idat[fl]; float qdat[fl];
     for(i=0; i<fl; i++) {
        c = pmt::to_complex(pmt::vector_ref(vector,i));
        idat[i] = c.real();
        qdat[i] = c.imag();
     }
     // 2nd tuple component is number of candidate frequencies
     int npk = pmt::to_long(pmt::tuple_ref(tuple,1));
     // 3rd tuple component is candidate frequencies
     pmt::pmt_t candidate_freqs(pmt::tuple_ref(tuple,2));
     for (i = 0; i < npk; i++) {
       pmt::pmt_t freq_tuple(pmt::vector_ref(candidate_freqs,i));
       candidates[i].m_type = (Modes)pmt::to_long(pmt::tuple_ref(freq_tuple,0));
       candidates[i].freq = pmt::to_double(pmt::tuple_ref(freq_tuple,1));
       candidates[i].snr = pmt::to_double(pmt::tuple_ref(freq_tuple,2));
       candidates[i].sync = pmt::to_double(pmt::tuple_ref(freq_tuple,3));
       candidates[i].shift = pmt::to_long(pmt::tuple_ref(freq_tuple,4));
       switch(candidates[i].m_type) {
         case linear : {
            candidates[i].m_linear.drift =
             pmt::to_double(pmt::tuple_ref(freq_tuple,5));
            break;
         }
         case nonlinear : {
            candidates[i].m_nonlinear.va =
             pmt::to_double(pmt::tuple_ref(freq_tuple,5)),
            candidates[i].m_nonlinear.ma =
             pmt::to_double(pmt::tuple_ref(freq_tuple,6)),
            candidates[i].m_nonlinear.xa =
             (int)pmt::to_long(pmt::tuple_ref(freq_tuple,7)),
            candidates[i].m_nonlinear.ya =
             (int)pmt::to_long(pmt::tuple_ref(freq_tuple,8));
            candidates[i].m_nonlinear.vb = 0;
            candidates[i].m_nonlinear.mb = 0;
            candidates[i].m_nonlinear.xb = 0;
            candidates[i].m_nonlinear.yb = 0;
            candidates[i].m_linear.drift = 0;
            break;
         }
       }
     }
     /*
      Refine the estimates of freq, shift using sync as a metric.
      Sync is calculated such that it is a float taking values in the range
      [0.0,1.0].
      Function sync_and_demodulate has three modes of operation
      mode is the last argument:
      0 = no frequency or drift search. find best time lag.
      1 = no time lag or drift search. find best frequency.
      2 = no frequency or time lag search. Calculate soft-decision
          symbols using passed frequency and shift.
     */
     for (j=0; j<npk; j++) {
       log(">>>  Demodulation");
       log("    Baseband freq is %2.2f Hz", candidates[j].freq);
       log("    (6 Hz) SNR is %2.2f dB", candidates[j].snr);
       if (candidates[j].m_type==linear) {
         log("    Linear drift is %2.2f Hz", candidates[j].m_linear.drift);
       } else {
         log("    Nonlinear drift  va:%2.2f ma:%2.2f xa:%d ya:%d",
           candidates[j].m_nonlinear.va,
           candidates[j].m_nonlinear.ma,
           candidates[j].m_nonlinear.xa,
           candidates[j].m_nonlinear.ya
         );
       }
       memset(symbols,0,sizeof(char)*nbits*2);
       f1=candidates[j].freq;
       drift1=candidates[j].m_linear.drift;
       shift1=candidates[j].shift;
       sync1=candidates[j].sync;
       // coarse-grid lag and freq search, then if sync>minsync1 continue
       fstep=0.0; ifmin=0; ifmax=0;
       lagmin=shift1-128;
       lagmax=shift1+128;
       lagstep=64;
       sync_and_demodulate(candidates[j], idat, qdat, npoints, symbols, &f1,
         ifmin, ifmax, fstep, &shift1,
         lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);
       fstep=0.25; ifmin=-2; ifmax=2;
       sync_and_demodulate(candidates[j], idat, qdat, npoints, symbols, &f1,
         ifmin, ifmax, fstep, &shift1,
         lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);
       // ---------------------------------------------------
       // In the linear case, refine frequency drift estimate
       // ---------------------------------------------------
       if (candidates[j].m_type==linear) {
         fstep=0.0; ifmin=0; ifmax=0;
         float driftp,driftm,syncp,syncm;
         driftp=drift1+0.5;
         sync_and_demodulate(candidates[j], idat, qdat, npoints, symbols, &f1,
           ifmin, ifmax, fstep, &shift1,
           lagmin, lagmax, lagstep, &driftp, symfac, &syncp, 1);
           driftm=drift1-0.5;
         sync_and_demodulate(candidates[j], idat, qdat, npoints, symbols, &f1,
           ifmin, ifmax, fstep, &shift1,
           lagmin, lagmax, lagstep, &driftm, symfac, &syncm, 1);
         if(syncp>sync1) {
           drift1=driftp;
           sync1=syncp;
         } else if (syncm>sync1) {
           drift1=driftm;
           sync1=syncm;
         }
       }
       // fine-grid lag and freq search
       if( sync1 > minsync1 ) {
         lagmin=shift1-32; lagmax=shift1+32; lagstep=16;
         sync_and_demodulate(candidates[j], idat, qdat, npoints, symbols, &f1,
           ifmin, ifmax, fstep, &shift1,
           lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);
         // fine search over frequency
         fstep=0.05; ifmin=-2; ifmax=2;
         sync_and_demodulate(candidates[j], idat, qdat, npoints, symbols, &f1,
           ifmin, ifmax, fstep, &shift1,
           lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);
         worth_a_try = 1;
       } else {
         worth_a_try = 0;
       }
       int idt=0, ii=0, jiggered_shift;
       float y,sq,rms;
       not_decoded=1;
       while ( worth_a_try && not_decoded && idt<=(128/iifac)) {
         ii=(idt+1)/2;
         if( idt%2 == 1 ) ii=-ii;
         ii=iifac*ii;
         jiggered_shift=shift1+ii;
         sync_and_demodulate(candidates[j], idat, qdat, npoints, symbols, &f1,
           ifmin, ifmax, fstep,
           &jiggered_shift, lagmin, lagmax, lagstep, &drift1, symfac,
           &sync1, 2);
         sq=0.0;
         for(i=0; i<162; i++) {
           y=(float)symbols[i] - 128.0;
           sq += y*y;
         }
         rms=sqrt(sq/162.0);
         if((sync1 > minsync2) && (rms > minrms)) {
           deinterleave(symbols);
           not_decoded = fano(&metric,&cycles,&maxnp,decdata,symbols,nbits,
             mettab,delta,maxcycles);
           //log("*** Result of Fano (not_decoded) is %d",not_decoded);
          }
          idt++;
       }
       if( worth_a_try && !not_decoded ) {
         for(i=0; i<11; i++) {
           if( decdata[i]>127 ) {
             message[i]=decdata[i]-256;
           } else {
             message[i]=decdata[i];
           }
         }
         // increment frame count
         framecount++;
         log("Frame: %d", framecount);
         log("    Baseband freq is %2.2f Hz", candidates[j].freq);
         log("    (6 Hz) SNR is %2.2f dB", candidates[j].snr);
         if (candidates[j].m_type==linear) {
           log("    Linear drift is %2.2f Hz", candidates[j].m_linear.drift);
         } else {
           log("    Nonlinear drift  va:%2.2f ma:%2.2f xa:%d ya:%d",
             candidates[j].m_nonlinear.va,
             candidates[j].m_nonlinear.ma,
             candidates[j].m_nonlinear.xa,
             candidates[j].m_nonlinear.ya
           );
         }
         // log message in a file
         printtime();
         fprintf(msglogfile, "Frame: %d\n", framecount);
         fprintf(msglogfile, "Baseband freq is %2.2f Hz\n", candidates[j].freq);
         fprintf(msglogfile, "(6 Hz) SNR is %2.2f dB\n", candidates[j].snr);
         if (candidates[j].m_type==linear) {
           fprintf(msglogfile,
             "Linear drift is %2.2f Hz\n", candidates[j].m_linear.drift);
         } else {
           fprintf(msglogfile,
             "Nonlinear drift  va:%2.2f ma:%2.2f xa:%d ya:%d\n",
             candidates[j].m_nonlinear.va,
             candidates[j].m_nonlinear.ma,
             candidates[j].m_nonlinear.xa,
             candidates[j].m_nonlinear.ya
           );
         }
         fprintf(msglogfile, "Data: ");
         logmessage(message); fprintf(msglogfile, "\n");
         fprintf(msglogfile, "\n");
         fflush (msglogfile);
         // push the packed message
         pmt::pmt_t payload(pmt::make_blob(&message[0], 7));
         pmt::pmt_t pdu(pmt::cons(pmt::PMT_NIL, payload));
         message_port_pub(out_port, pdu);
       }
     } // end of for loop
     return;
    }
   }   /* namespace uwspr */
} /* namespace gr */
