/* -*- c++ -*- */
/*
 * Copyright 2017 Michel Barbeau, Carleton University.
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
#include "decoder_impl.h"
#include <stdio.h>
#include <stdlib.h>
#include <boost/thread/thread.hpp>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <stdbool.h>
#include <volk/volk.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
//#include "wsprd.h"
#include "fano.h"

namespace gr {
namespace uwspr {

decoder::sptr
decoder::make(int maxdrift)
{
   return gnuradio::get_initial_sptr
             (new decoder_impl(maxdrift));
}

/*
 * The private constructor
 */
decoder_impl::decoder_impl(int maxdrift)
   : gr::sync_block("decoder",
                    gr::io_signature::make(1, 1, sizeof(gr_complex)),
                    gr::io_signature::make(0, 0, 0))
{
   // maximum frequency drift
   this->maxdrift = maxdrift;
   // create and register the asynchronous output port
   out_port = pmt::mp("out");
   message_port_register_out(out_port);
   // constants
   fs = 375;      // sampling frequency (samples/second)
   fl = 45000;
   count = 0;
   const int alignment_multiple =
      volk_get_alignment() / sizeof(float);
   set_alignment(std::max(1,alignment_multiple));
   buffer.set_capacity(2*fl);
   sem_init(&semaphore, 0, 0);
   // create decoding thread(s)
   boost::thread * t;
   for (int i=0; i<10; ++i) {
      t = new boost::thread(boost::bind(&decoder_impl::decode, this));
      t->detach();
   }
   // WSPR initialization
   strcpy(wisdom_fname,".");
   if(data_dir != NULL) {
      strcpy(wisdom_fname,data_dir);
   }
   strncat(wisdom_fname,"/wspr_wisdom.dat",20);
   if ((fp_fftwf_wisdom_file = fopen(wisdom_fname, "r"))) { //Open FFTW wisdom
      fftwf_import_wisdom_from_file(fp_fftwf_wisdom_file);
      fclose(fp_fftwf_wisdom_file);
   }
   tfano=0.0; treadwav=0.0; tcandidates=0.0; tsync0=0.0;
   tsync1=0.0; tsync2=0.0; ttotal=0.0;
   // stack decoding initialization
   stacksize=200000;
   stackdecoder=0;
   if(stackdecoder) {
      stack=(struct snode *)malloc(stacksize*sizeof(struct snode));
   }
   // metric table initialization
#include "./metric_tables.c"
   bias=0.45;      //Fano metric bias (used for both Fano and stack algorithms)
   for(int i=0; i<256; i++) {
      mettab[0][i]=round( 10*(metric_tables[2][i]-bias) );
      mettab[1][i]=round( 10*(metric_tables[2][255-i]-bias) );
   }
   // initialization of date and time
   date[0]='\0';
   uttime[0]='\0';
   // --- preliminaries for one-dimensional DFT of size 512
   // Do windowed ffts over 2 symbols, stepped by half symbols
   // Note: One symbol (baud) is 256 samples. 2 symbols times
   // 256 samples/symbol is 512 samples. "npoints" divided by
   // 512 is the number 2 symbol chunks. 4 times the number of symbols
   // is the number of half symbol steps (i.e., nffts).
   npoints = 45000;
   nffts=4*floor(npoints/512)-1;
   // allocate memory for FFT input array
   fftin=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);
   // allocate memory for FFT output array
   fftout=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*512);
   // plan allocation (create object containing data needed by FFTW)
   // 512 = size of transform
   // fftin = pointer to input array
   // fftout = pointer to output array
   // FFTW_FORWARD (-1) = direction of transform (sign of exponent in transform)
   // PATIENCE = flags
   PLAN3 = fftwf_plan_dft_1d(512, fftin, fftout, FFTW_FORWARD, PATIENCE);
   // generate a windowing function (half-sine wave cycle [0 pi]), in 512 steps
   for(int i=0; i<512; i++) {
      w[i]=sin(0.006147931*i);
   }
}

/*
 * Our virtual destructor.
 */
decoder_impl::~decoder_impl()
{
   fftwf_destroy_plan(PLAN3);
   if( stackdecoder ) {
      free(stack);
   }
   fftwf_free(fftin);
   fftwf_free(fftout);
}

// WSPR Code

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
 *
 * Storing only one symbol per byte uses more space, but it is faster
 * and easier than trying to pack them more compactly.
 */
int decoder_impl::encode(
   unsigned char *symbols,      // Output buffer, 2*nbytes*8
   unsigned char *data,      // Input buffer, nbytes
   unsigned int nbytes)      // Number of bytes in data
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
 *
 * Decode packet with the Fano algorithm.
 * Return 0 on success, -1 on timeout
 *
 * Decode packet with the Fano algorithm.
 * Return 0 on success, -1 on timeout
 */
int decoder_impl::fano(
   unsigned int  *metric,      // Final path metric (returned value)
   unsigned int  *cycles,      // Cycle count (returned value)
   unsigned int  *maxnp,      // Progress before timeout (returned value)
   unsigned char *data,      // Decoded output data
   unsigned char *symbols,      // Raw deinterleaved input symbols
   unsigned int nbits,      // Number of output bits
   int mettab[2][256],      // Metric table, [sent sym][rx symbol]
   int delta,      // Threshold adjust parameter
   unsigned int maxcycles)      // Decoding timeout in cycles per bit
{
   struct node *nodes;      // First node
   struct node *np;            // Current node
   struct node *lastnode;      // Last node
   struct node *tail;      // First node of tail
   int t;        // Threshold
   int m0,m1;
   int ngamma;
   unsigned int lsym;
   unsigned int i;

   if((nodes = (struct node *)malloc((nbits+1)*sizeof(struct node))) == NULL) {
      printf("malloc failed\n");
      return 0;
   }
   lastnode = &nodes[nbits-1];
   tail = &nodes[nbits-31];
   *maxnp = 0;

/* Compute all possible branch metrics for each symbol pair
 * This is the only place we actually look at the raw input symbols
 */
   for(np=nodes; np <= lastnode; np++) {
      np->metrics[0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
      np->metrics[1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
      np->metrics[2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
      np->metrics[3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
      symbols += 2;
   }
   np = nodes;
   np->encstate = 0;

// Compute and sort branch metrics from root node */
   ENCODE(lsym,np->encstate);      // 0-branch (LSB is 0)
   m0 = np->metrics[lsym];

/* Now do the 1-branch. To save another ENCODE call here and
 * inside the loop, we assume that both polynomials are odd,
 * providing complementary pairs of branch symbols.

 * This code should be modified if a systematic code were used.
 */

   m1 = np->metrics[3^lsym];
   if(m0 > m1) {
      np->tm[0] = m0;                           // 0-branch has better metric
      np->tm[1] = m1;
   } else {
      np->tm[0] = m1;                           // 1-branch is better
      np->tm[1] = m0;
      np->encstate++;                       // Set low bit
   }
   np->i = 0;                                 // Start with best branch
   maxcycles *= nbits;
   np->gamma = t = 0;

   // Start the Fano decoder
   for(i=1; i <= maxcycles; i++) {
      if((int)(np-nodes) > (int)*maxnp) *maxnp=(int)(np-nodes);
#ifdef  debug
      printf("k=%ld, g=%ld, t=%d, m[%d]=%d, maxnp=%d, encstate=%lx\n",
             np-nodes,np->gamma,t,np->i,np->tm[np->i],*maxnp,np->encstate);
#endif
// Look forward */
      ngamma = np->gamma + np->tm[np->i];
      if(ngamma >= t) {
         if(np->gamma < t + delta) {                // Node is acceptable
            /* First time we've visited this node;
             * Tighten threshold.
             *
             * This loop could be replaced with
             *   t += delta * ((ngamma - t)/delta);
             * but the multiply and divide are slower.
             */
            while(ngamma >= t + delta) t += delta;
         }
         np[1].gamma = ngamma;                  // Move forward
         np[1].encstate = np->encstate << 1;
         if( ++np == (lastnode+1) ) {
            break;                        // Done!
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
               np->tm[0] = m0;                          // 0-branch is better
               np->tm[1] = m1;
            } else {
               np->tm[0] = m1;                          // 1-branch is better
               np->tm[1] = m0;
               np->encstate++;                          // Set low bit
            }
         }
         np->i = 0;                       // Start with best branch
         continue;
      }
      // Threshold violated, can't go forward
      for(;; ) {                                // Look backward
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
            np->i++;                     // Search next best branch
            np->encstate ^= 1;
            break;
         }                                // else keep looking back
      }
   }
   *metric =  np->gamma;                  // Return the final path metric

   // Copy decoded data to user's buffer
   nbits >>= 3;
   np = &nodes[7];
   while(nbits-- != 0) {
      *data++ = np->encstate;
      np += 8;
   }
   *cycles = i+1;

   free(nodes);
   if(i >= maxcycles) return -1;          // Decoder timed out
   return 0;                    // Successful completion
}

// synchronization bit pattern
unsigned char pr3[162]=
{1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
 0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
 0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
 1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
 0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
 0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
 0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
 0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
 0,0};

void decoder_impl::sync_and_demodulate(
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
//
// Input signal:
//    id = in-phase (float array of maxpts=65536)
//    qd = quarature (float array of maxpts=65536)
//    np = number of points
// symbols = (char array of size nbits (81) times 2)
   static float fplast=-10000.0;
   // delta time (dt) and delta (frequency)
   static float dt=1.0/375.0, df=375.0/256.0;
   static float pi=3.14159265358979323846;
   float twopidt, df15=df*1.5, df05=df*0.5;

   int i, j, k, lag;
   float i0[162],q0[162],i1[162],q1[162],i2[162],q2[162],i3[162],q3[162];
   float p0,p1,p2,p3,cmet,totp,syncmax,fac;
   float c0[256],s0[256],c1[256],s1[256],c2[256],s2[256],c3[256],s3[256];
   float dphi0, cdphi0, sdphi0, dphi1, cdphi1, sdphi1, dphi2, cdphi2, sdphi2,
         dphi3, cdphi3, sdphi3;
   float ss; // correlation with synchronization soft symbols
   float f0=0.0, fp, fbest=0.0, fsum=0.0, f2sum=0.0, fsymb[162];
   int best_shift = 0, ifreq;

   syncmax=-1e30;
   if( mode == 0 ) {ifmin=0; ifmax=0; fstep=0.0; f0=*f1; }
   if( mode == 1 ) {lagmin=*shift1; lagmax=*shift1; f0=*f1; }
   if( mode == 2 ) {lagmin=*shift1; lagmax=*shift1; ifmin=0; ifmax=0; f0=*f1; }

   twopidt=2*pi*dt;
   for(ifreq=ifmin; ifreq<=ifmax; ifreq++) { // frequency search
      f0=*f1+ifreq*fstep;
      for(lag=lagmin; lag<=lagmax; lag=lag+lagstep) { // time search
         ss=0.0;
         totp=0.0;
         for (i=0; i<162; i++) {
            fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;  // drift search
            // do for each symbol
            if( i==0 || (fp != fplast) ) { // only calculate sin/cos if necessary
               dphi0=twopidt*(fp-df15); // phase of symbol 0
               cdphi0=cos(dphi0); // cosine increment of a sample
               sdphi0=sin(dphi0); // sine increment of a sample

               dphi1=twopidt*(fp-df05); // phase of symbol 1
               cdphi1=cos(dphi1);
               sdphi1=sin(dphi1);

               dphi2=twopidt*(fp+df05); // phase of symbol 2
               cdphi2=cos(dphi2);
               sdphi2=sin(dphi2);

               dphi3=twopidt*(fp+df15); // phase of symbol 3
               cdphi3=cos(dphi3);
               sdphi3=sin(dphi3);

               c0[0]=1; s0[0]=0; // start at phase 0 radian
               c1[0]=1; s1[0]=0;
               c2[0]=1; s2[0]=0;
               c3[0]=1; s3[0]=0;

               // calculate cosine and sine for each sample of symbol
               for (j=1; j<256; j++) {
                  c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                  s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                  c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                  s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                  c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                  s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                  c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                  s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
               }
               fplast = fp;
            }

            i0[i]=0.0; q0[i]=0.0;
            i1[i]=0.0; q1[i]=0.0;
            i2[i]=0.0; q2[i]=0.0;
            i3[i]=0.0; q3[i]=0.0;

            // correlation
            for (j=0; j<256; j++) {
               k=lag+i*256+j;
               if( (k>0) && (k<np) ) {
                  i0[i]=i0[i] + id[k]*c0[j] + qd[k]*s0[j];
                  q0[i]=q0[i] - id[k]*s0[j] + qd[k]*c0[j];
                  i1[i]=i1[i] + id[k]*c1[j] + qd[k]*s1[j];
                  q1[i]=q1[i] - id[k]*s1[j] + qd[k]*c1[j];
                  i2[i]=i2[i] + id[k]*c2[j] + qd[k]*s2[j];
                  q2[i]=q2[i] - id[k]*s2[j] + qd[k]*c2[j];
                  i3[i]=i3[i] + id[k]*c3[j] + qd[k]*s3[j];
                  q3[i]=q3[i] - id[k]*s3[j] + qd[k]*c3[j];
               }
            }
            p0=i0[i]*i0[i] + q0[i]*q0[i];
            p1=i1[i]*i1[i] + q1[i]*q1[i];
            p2=i2[i]*i2[i] + q2[i]*q2[i];
            p3=i3[i]*i3[i] + q3[i]*q3[i];

            p0=sqrt(p0);
            p1=sqrt(p1);
            p2=sqrt(p2);
            p3=sqrt(p3);

            totp=totp+p0+p1+p2+p3; // total power
            cmet=(p1+p3)-(p0+p2); // metric
            ss = (pr3[i] == 1) ? ss+cmet : ss-cmet;
            if( mode == 2) { //Compute soft symbols (ss)
               if(pr3[i]==1) {
                  // synchro bit is 1
                  fsymb[i]=p3-p1; // difference of powers at data bits 1 and 0
               } else {
                  // synchro bit is 0
                  fsymb[i]=p2-p0; // difference of powers at data bits 1 and 0
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

/******************************************************************************
   Fully coherent signal subtraction
 *******************************************************************************/
void decoder_impl::subtract_signal2(float *id, float *qd, long np,
                                    float f0, int shift0, float drift0, unsigned char* channel_symbols)
{
   float dt=1.0/375.0, df=375.0/256.0;
   float pi=4.*atan(1.0), twopidt, phi=0, dphi, cs;
   int i, j, k, ii, nsym=162, nspersym=256,  nfilt=256;      //nfilt must be even number.
   int nsig=nsym*nspersym;
   int nc2=45000;

   float *refi, *refq, *ci, *cq, *cfi, *cfq;

   refi=(float *)malloc(sizeof(float)*nc2);
   refq=(float *)malloc(sizeof(float)*nc2);
   ci=(float *)malloc(sizeof(float)*nc2);
   cq=(float *)malloc(sizeof(float)*nc2);
   cfi=(float *)malloc(sizeof(float)*nc2);
   cfq=(float *)malloc(sizeof(float)*nc2);

   memset(refi,0,sizeof(float)*nc2);
   memset(refq,0,sizeof(float)*nc2);
   memset(ci,0,sizeof(float)*nc2);
   memset(cq,0,sizeof(float)*nc2);
   memset(cfi,0,sizeof(float)*nc2);
   memset(cfq,0,sizeof(float)*nc2);

   twopidt=2.0*pi*dt;

   /******************************************************************************
      Measured signal:                    s(t)=a(t)*exp( j*theta(t) )
      Reference is:                       r(t) = exp( j*phi(t) )
      Complex amplitude is estimated as:  c(t)=LPF[s(t)*conjugate(r(t))]
      so c(t) has phase angle theta-phi
      Multiply r(t) by c(t) and subtract from s(t), i.e. s'(t)=s(t)-c(t)r(t)
    *******************************************************************************/

   // create reference wspr signal vector, centered on f0.
   //
   for (i=0; i<nsym; i++) {

      cs=(float)channel_symbols[i];

      dphi=twopidt*
            (
         f0 + (drift0/2.0)*((float)i-(float)nsym/2.0)/((float)nsym/2.0)
         + (cs-1.5)*df
            );

      for ( j=0; j<nspersym; j++ ) {
         ii=nspersym*i+j;
         refi[ii]=cos(phi);                //cannot precompute sin/cos because dphi is changing
         refq[ii]=sin(phi);
         phi=phi+dphi;
      }
   }

   // s(t) * conjugate(r(t))
   // beginning of first symbol in reference signal is at i=0
   // beginning of first symbol in received data is at shift0.
   // filter transient lasts nfilt samples
   // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
   for (i=0; i<nsym*nspersym; i++) {
      k=shift0+i;
      if( (k>0) && (k<np) ) {
         ci[i+nfilt] = id[k]*refi[i] + qd[k]*refq[i];
         cq[i+nfilt] = qd[k]*refi[i] - id[k]*refq[i];
      }
   }

   //lowpass filter and remove startup transient
   float w[nfilt], norm=0, partialsum[nfilt];
   memset(partialsum,0,sizeof(float)*nfilt);
   for (i=0; i<nfilt; i++) {
      w[i]=sin(pi*(float)i/(float)(nfilt-1));
      norm=norm+w[i];
   }
   for (i=0; i<nfilt; i++) {
      w[i]=w[i]/norm;
   }
   for (i=1; i<nfilt; i++) {
      partialsum[i]=partialsum[i-1]+w[i];
   }

   // LPF
   for (i=nfilt/2; i<45000-nfilt/2; i++) {
      cfi[i]=0.0; cfq[i]=0.0;
      for (j=0; j<nfilt; j++) {
         cfi[i]=cfi[i]+w[j]*ci[i-nfilt/2+j];
         cfq[i]=cfq[i]+w[j]*cq[i-nfilt/2+j];
      }
   }

   // subtract c(t)*r(t) here
   // (ci+j*cq)(refi+j*refq)=(ci*refi-cq*refq)+j(ci*refq)+cq*refi)
   // beginning of first symbol in reference signal is at i=nfilt
   // beginning of first symbol in received data is at shift0.
   for (i=0; i<nsig; i++) {
      if( i<nfilt/2 ) {           // take care of the end effect (LPF step response) here
         norm=partialsum[nfilt/2+i];
      } else if( i>(nsig-1-nfilt/2) ) {
         norm=partialsum[nfilt/2+nsig-1-i];
      } else {
         norm=1.0;
      }
      k=shift0+i;
      j=i+nfilt;
      if( (k>0) && (k<np) ) {
         id[k]=id[k] - (cfi[j]*refi[i]-cfq[j]*refq[i])/norm;
         qd[k]=qd[k] - (cfi[j]*refq[i]+cfq[j]*refi[i])/norm;
      }
   }

   free(refi);
   free(refq);
   free(ci);
   free(cq);
   free(cfi);
   free(cfq);

   return;
}

unsigned long decoder_impl::writec2file(char *c2filename, int trmin,
                                        double freq, float *idat, float *qdat)
{
   int i;
   float *buffer;
   buffer=(float *)malloc(sizeof(float)*2*45000);
   memset(buffer,0,sizeof(float)*2*45000);

   FILE *fp;

   fp = fopen(c2filename,"wb");
   if( fp == NULL ) {
      fprintf(stderr, "Could not open c2 file '%s'\n", c2filename);
      free(buffer);
      return 0;
   }
   unsigned long nwrite = fwrite(c2filename,sizeof(char),14,fp);
   nwrite = fwrite(&trmin, sizeof(int), 1, fp);
   nwrite = fwrite(&freq, sizeof(double), 1, fp);

   for(i=0; i<45000; i++) {
      buffer[2*i]=idat[i];
      buffer[2*i+1]=-qdat[i];
   }

   nwrite = fwrite(buffer, sizeof(float), 2*45000, fp);
   if( nwrite == 2*45000 ) {
      return nwrite;
   } else {
      free(buffer);
      return 0;
   }
}

/* WSPR uses the Layland-Lushbaugh code
 * Nonsystematic, non-quick look-in, dmin=?, dfree=?
 */
#define POLY1 0xf2d05351
#define POLY2 0xe4613c47

//Decoder - returns 0 on success, -1 on timeout
int decoder_impl::jelinek(
   unsigned int *metric,          /* Final path metric (returned value) */
   unsigned int *cycles,          /* Cycle count (returned value) */
   unsigned char *data,           /* Decoded output data */
   unsigned char *symbols,          /* Raw deinterleaved input symbols */
   unsigned int nbits,          /* Number of output bits */
   unsigned int stacksize,
   struct snode *stack,
   int mettab[2][256],          /* Metric table, [sent sym][rx symbol] */
   unsigned int maxcycles)         /* Decoding timeout in cycles per bit */
{

   // Compute branch metrics for each symbol pair
   // The sequential decoding algorithm only uses the metrics, not the
   // symbol values.
   unsigned int i;
   long int metrics[81][4];
   for(i=0; i<nbits; i++) {
      metrics[i][0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
      metrics[i][1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
      metrics[i][2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
      metrics[i][3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
      symbols += 2;
   }

   // zero the stack
   memset(stack,0,stacksize*sizeof(struct snode));

   // initialize the loop variables
   unsigned int lsym, ntail=31;
   uint64_t encstate=0;
   unsigned int nbuckets=1000;
   unsigned int low_bucket=nbuckets-1;      //will be set on first run-through
   unsigned int high_bucket=0;
   unsigned int *buckets, bucket;
   buckets=(unsigned int *)malloc(nbuckets*sizeof(unsigned int));
   memset(buckets,0,nbuckets*sizeof(unsigned int));
   unsigned int ptr=1;
   unsigned int stackptr=1;      //pointer values of 0 are reserved (they mean that a bucket is empty)
   unsigned int depth=0, nbits_minus_ntail=nbits-ntail;
   unsigned int stacksize_minus_1=stacksize-1;
   long int totmet0, totmet1, gamma=0;

   unsigned int ncycles=maxcycles*nbits;
   /********************* Start the stack decoder *****************/
   for (i=1; i <= ncycles; i++) {
#ifdef DEBUG
      printf("***stackptr=%ld, depth=%d, gamma=%d, encstate=%lx, bucket %d, low_bucket %d, high_bucket %d\n",
             stackptr, depth, gamma, encstate, bucket, low_bucket, high_bucket);
#endif
      // no need to store more than 7 bytes (56 bits) for encoder state because
      // only 50 bits are not 0's.
      if( depth < 56 ) {
         encstate=encstate<<1;
         ENCODE(lsym,encstate);                // get channel symbols associated with the 0 branch
      } else {
         ENCODE(lsym,encstate<<(depth-55));
      }

      // lsym are the 0-branch channel symbols and 3^lsym are the 1-branch
      // channel symbols (due to a special property of our generator polynomials)
      totmet0 = gamma+metrics[depth][lsym];           // total metric for 0-branch daughter node
      totmet1 = gamma+metrics[depth][3^lsym];           // total metric for 1-branch daughter node
      depth++;           //the depth of the daughter nodes

      bucket=(totmet0>>5)+200;           //fast, but not particularly safe - totmet can be negative
      if( bucket > high_bucket ) high_bucket=bucket;
      if( bucket < low_bucket ) low_bucket=bucket;

      // place the 0 node on the stack, overwriting the parent (current) node
      stack[ptr].encstate=encstate;
      stack[ptr].gamma=totmet0;
      stack[ptr].depth=depth;
      stack[ptr].jpointer=buckets[bucket];
      buckets[bucket]=ptr;

      // if in the tail, only need to evaluate the "0" branch.
      // Otherwise, enter this "if" and place the 1 node on the stack,
      if( depth <= nbits_minus_ntail ) {
         if( stackptr < stacksize_minus_1 ) {
            stackptr++;
            ptr=stackptr;
         } else {                // stack full
            while( buckets[low_bucket] == 0 ) {                     //write latest to where the top of the lowest bucket points
               low_bucket++;
            }
            ptr=buckets[low_bucket];
            buckets[low_bucket]=stack[ptr].jpointer;                     //make bucket point to next older entry
         }

         bucket=(totmet1>>5)+200;                //this may not be safe on all compilers
         if( bucket > high_bucket ) high_bucket=bucket;
         if( bucket < low_bucket ) low_bucket=bucket;

         stack[ptr].encstate=encstate+1;
         stack[ptr].gamma=totmet1;
         stack[ptr].depth=depth;
         stack[ptr].jpointer=buckets[bucket];
         buckets[bucket]=ptr;
      }

      // pick off the latest entry from the high bucket
      while( buckets[high_bucket] == 0 ) {
         high_bucket--;
      }
      ptr=buckets[high_bucket];
      buckets[high_bucket]=stack[ptr].jpointer;
      depth=stack[ptr].depth;
      gamma=stack[ptr].gamma;
      encstate=stack[ptr].encstate;

      // we are done if the top entry on the stack is at depth nbits
      if (depth == nbits) {
         break;
      }
   }

   *cycles = i+1;
   *metric =  gamma;      /* Return final path metric */

   //    printf("cycles %d stackptr=%d, depth=%d, gamma=%d, encstate=%lx\n",
   //           *cycles, stackptr, depth, *metric, encstate);

   for (i=0; i<7; i++) {
      data[i]=(encstate>>(48-i*8))&(0x00000000000000ff);
   }
   for (i=7; i<11; i++) {
      data[i]=0;
   }

   if(*cycles/nbits >= maxcycles)      //timed out
   {
      return -1;
   }
   return 0;      //success
}

int decoder_impl::floatcomp(const void* elem1, const void* elem2)
{
   if(*(const float*)elem1 < *(const float*)elem2)
      return -1;
   return *(const float*)elem1 > *(const float*)elem2;
}

void decoder_impl::deinterleave(unsigned char *sym)
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

char decoder_impl::get_locator_character_code(char ch) {
   if( ch >=48 && ch <=57 ) {      //0-9
      return ch-48;
   }
   if( ch == 32 ) {      //space
      return 36;
   }
   if( ch >= 65 && ch <= 82 ) {      //A-Z
      return ch-65;
   }
   return -1;
}

char decoder_impl::get_callsign_character_code(char ch) {
   if( ch >=48 && ch <=57 ) {      //0-9
      return ch-48;
   }
   if( ch == 32 ) {      //space
      return 36;
   }
   if( ch >= 65 && ch <= 90 ) {      //A-Z
      return ch-55;
   }
   return -1;
}

long unsigned int decoder_impl::pack_grid4_power(char *grid4, int power) {
   long unsigned int m;

   m=(179-10*grid4[0]-grid4[2])*180+10*grid4[1]+grid4[3];
   m=m*128+power+64;
   return m;
}

long unsigned int decoder_impl::pack_call(char *callsign) {
   unsigned int i;
   long unsigned int n;
   char call6[6];
   memset(call6,32,sizeof(char)*6);
   // callsign is 6 characters in length. Exactly.
   size_t call_len = strlen(callsign);
   if( call_len > 6 ) {
      return 0;
   }
   if( isdigit(*(callsign+2)) ) {
      for (i=0; i<call_len; i++) {
         if( callsign[i] == 0 ) {
            call6[i]=32;
         } else {
            call6[i]=*(callsign+i);
         }
      }
   } else if( isdigit(*(callsign+1)) ) {
      call6[0]=32;
      for (i=1; i<call_len+1; i++) {
         if( callsign[i-1]==0 ) {
            call6[i]=32;
         } else {
            call6[i]=*(callsign+i-1);
         }
      }
   }
   for (i=0; i<6; i++) {
      call6[i]=get_callsign_character_code(call6[i]);
   }
   n = call6[0];
   n = n*36+call6[1];
   n = n*10+call6[2];
   n = n*27+call6[3]-10;
   n = n*27+call6[4]-10;
   n = n*27+call6[5]-10;
   return n;
}

void decoder_impl::pack_prefix(char *callsign, int32_t *n, int32_t *m, int32_t *nadd ) {
   size_t i;
   char *call6;
   call6=(char *)malloc(sizeof(char)*6);
   memset(call6,32,sizeof(char)*6);
   size_t i1=strcspn(callsign,"/");

   if( callsign[i1+2] == 0 ) {
      //single char suffix
      for (i=0; i<i1; i++) {
         call6[i]=callsign[i];
      }
      *n=pack_call(call6);
      *nadd=1;
      int nc = callsign[i1+1];
      if( nc >= 48 && nc <= 57 ) {
         *m=nc-48;
      } else if ( nc >= 65 && nc <= 90 ) {
         *m=nc-65+10;
      } else {
         *m=38;
      }
      *m=60000-32768+*m;
   } else if( callsign[i1+3]==0 ) {
      //two char suffix
      for (i=0; i<i1; i++) {
         call6[i]=callsign[i];
      }
      *n=pack_call(call6);
      *nadd=1;
      *m=10*(callsign[i1+1]-48)+(callsign[i1+2]-48);
      *m=60000 + 26 + *m;
   } else {
      char* pfx=strtok(callsign,"/");
      call6=strtok(NULL," ");
      *n=pack_call(call6);
      size_t plen=strlen(pfx);
      if( plen ==1 ) {
         *m=36;
         *m=37*(*m)+36;
      } else if( plen == 2 ) {
         *m=36;
      } else {
         *m=0;
      }
      for (i=0; i<plen; i++) {
         int nc = callsign[i];
         if( nc >= 48 && nc <= 57 ) {
            nc=nc-48;
         } else if ( nc >= 65 && nc <= 90 ) {
            nc=nc-65+10;
         } else {
            nc=36;
         }
         *m=37*(*m)+nc;
      }
      *nadd=0;
      if( *m > 32768 ) {
         *m=*m-32768;
         *nadd=1;
      }
   }
}

void decoder_impl::interleave(unsigned char *sym)
{
   unsigned char tmp[162];
   unsigned char p, i, j;

   p=0;
   i=0;
   while (p<162) {
      j=((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
      if (j < 162 ) {
         tmp[j]=sym[p];
         p=p+1;
      }
      i=i+1;
   }
   for (i=0; i<162; i++) {
      sym[i]=tmp[i];
   }
}

int decoder_impl::get_wspr_channel_symbols(char* rawmessage, char* hashtab, unsigned char* symbols) {
   int m=0, n=0, ntype=0;
   int i, j, ihash;
   unsigned char pr3[162]=
   {1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
    0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
    0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
    1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
    0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
    0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
    0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
    0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
    0,0};
   int nu[10]={0,-1,1,0,-1,2,1,0,-1,1};
   char *callsign, *grid, *powstr;
   char grid4[5], message[23];

   memset(message,0,sizeof(char)*23);
   i=0;
   while ( rawmessage[i] != 0 && i<23 ) {
      message[i]=rawmessage[i];
      i++;
   }

   size_t i1=strcspn(message," ");
   size_t i2=strcspn(message,"/");
   size_t i3=strcspn(message,"<");
   size_t i4=strcspn(message,">");
   size_t mlen=strlen(message);

   // Use the presence and/or absence of "<" and "/" to decide what
   // type of message. No sanity checks! Beware!

   if( (i1>3) & (i1<7) & (i2==mlen) & (i3==mlen) ) {
      // Type 1 message: K9AN EN50 33
      //                 xxnxxxx xxnn nn
      callsign = strtok(message," ");
      grid = strtok(NULL," ");
      powstr = strtok(NULL," ");
      int power = atoi(powstr);
      n = pack_call(callsign);

      for (i=0; i<4; i++) {
         grid4[i]=get_locator_character_code(*(grid+i));
      }
      m = pack_grid4_power(grid4,power);

   } else if ( i3 == 0 && i4 < mlen ) {
      // Type 3:      <K1ABC> EN50WC 33
      //          <PJ4/K1ABC> FK52UD 37
      // send hash instead of callsign to make room for 6 char grid.
      // if 4-digit locator is specified, 2 spaces are added to the end.
      callsign=strtok(message,"<> ");
      grid=strtok(NULL," ");
      powstr=strtok(NULL," ");
      int power = atoi(powstr);
      if( power < 0 ) power=0;
      if( power > 60 ) power=60;
      power=power+nu[power%10];
      ntype=-(power+1);
      ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
      m=128*ihash + ntype + 64;

      char grid6[6];
      memset(grid6,32,sizeof(char)*6);
      j=strlen(grid);
      for(i=0; i<j-1; i++) {
         grid6[i]=grid[i+1];
      }
      grid6[5]=grid[0];
      n=pack_call(grid6);
   } else if ( i2 < mlen ) {      // just looks for a right slash
      // Type 2: PJ4/K1ABC 37
      callsign=strtok(message," ");
      if( strlen(callsign) < i2 ) return 0;            //guards against pathological case
      powstr=strtok(NULL," ");
      int power = atoi(powstr);
      if( power < 0 ) power=0;
      if( power > 60 ) power=60;
      power=power+nu[power%10];
      int n1, ng, nadd;
      pack_prefix(callsign, &n1, &ng, &nadd);
      ntype=power + 1 + nadd;
      m=128*ng+ntype+64;
      n=n1;
   } else {
      return 0;
   }

   // pack 50 bits + 31 (0) tail bits into 11 bytes
   unsigned char it, data[11];
   memset(data,0,sizeof(char)*11);
   it=0xFF & (n>>20);
   data[0]=it;
   it=0xFF & (n>>12);
   data[1]=it;
   it=0xFF & (n>>4);
   data[2]=it;
   it= ((n&(0x0F))<<4) + ((m>>18)&(0x0F));
   data[3]=it;
   it=0xFF & (m>>10);
   data[4]=it;
   it=0xFF & (m>>2);
   data[5]=it;
   it=(m & 0x03)<<6;
   data[6]=it;
   data[7]=0;
   data[8]=0;
   data[9]=0;
   data[10]=0;

   //if( printdata ) {
   //    printf("Data is :");
   //    for (i=0; i<11; i++) {
   //        printf("%02X ",data[i]);
   //    }
   //    printf("\n");
   //}

   // make sure that the 11-byte data vector is unpackable
   // unpack it with the routine that the decoder will use and display
   // the result. let the operator decide whether it worked.

   char *check_call_loc_pow, *check_callsign;
   check_call_loc_pow=(char *)malloc(sizeof(char)*23);
   check_callsign=(char *)malloc(sizeof(char)*13);
   signed char check_data[11];
   memcpy(check_data,data,sizeof(char)*11);
   //unpk_(check_data,hashtab,check_call_loc_pow,check_callsign);
   //    printf("Will decode as: %s\n",check_call_loc_pow);

   unsigned int nbytes=11;      // The message with tail is packed into almost 11 bytes.
   unsigned char channelbits[nbytes*8*2];      /* 162 rounded up */
   memset(channelbits,0,sizeof channelbits);

   encode(channelbits,data,nbytes);

   interleave(channelbits);

   for (i=0; i<162; i++) {
      symbols[i]=2*channelbits[i]+pr3[i];
   }

   return 1;
}

#define max(x,y) ((x) > (y) ? (x) : (y))

// WSPR decoder
// idat = in-phase of input data
// qdat = quadrature of input data
// decodes = array of decoded messages
int decoder_impl::wsprd(
  float * idat, float * qdat,  struct result decodes[], int *uniques_ptr) {
   int i,j,k;
   int uniques=0;
   unsigned char *symbols, *decdata, *channel_symbols;
   signed char message[]={-9,13,-35,123,57,-39,64,0,0,0,0};
   char *callsign, *call_loc_pow;
   int c,delta,maxpts=65536,verbose=0,quickmode=0,more_candidates=0;
   int writenoise=0,wspr_type=2, ipass;
   int writec2=0, npasses=2, subtraction=1;
   int shift1, lagmin, lagmax, lagstep, ifmin, ifmax, worth_a_try, not_decoded;
   unsigned int nbits=81;
   unsigned int metric, cycles, maxnp;
   float df=375.0/256.0/2; // or df=775/256
   float freq0[200],snr0[200],drift0[200],sync0[200];
   int shift0[200];
   float dt=1.0/375.0, dt_print;
   double dialfreq_cmdline=0.0, dialfreq, freq_print;
   double dialfreq_error=0.0;
   float fmin=-110, fmax=110;
   float f1, fstep, sync1, drift1;
   float psavg[512];
   clock_t t0,t00;
   symbols=(unsigned char *)malloc(sizeof(char)*nbits*2);
   decdata=(unsigned char *)malloc(sizeof(char)*11);
   channel_symbols=(unsigned char *)malloc(sizeof(char)*nbits*2);
   callsign=(char *)malloc(sizeof(char)*13);
   call_loc_pow=(char *)malloc(sizeof(char)*23);
   float allfreqs[100];
   char allcalls[100][13];
   memset(allfreqs,0,sizeof(float)*100);
   memset(allcalls,0,sizeof(char)*100*13);
   int noprint=0, ndecodes_pass=0;
   // Parameters used for performance-tuning:
   unsigned int maxcycles=10000;               //Decoder timeout limit
   float minsync1=0.10;                        //First sync limit
   float minsync2=0.12;                        //Second sync limit
   int iifac=8;                                //Step size in final DT peakup
   int symfac=50;                              //Soft-symbol normalizing factor
   float minrms=52.0 * (symfac/64.0);         //Final test for plausible decoding
   delta=60;                                   //Fano threshold step
   t00=clock();
   // power spectrum representation
   float ps[512][nffts];      // power spectrum
   // timing
   struct timeval start, end;
   gettimeofday(&start, NULL);

//*************** main loop starts here *****************
   for (ipass=0; ipass<npasses; ipass++) {
      // printf("Pass index : %d\n",ipass);
      // one pass
      if( ipass > 0 && ndecodes_pass == 0 ) break;
      ndecodes_pass=0;
      memset(ps,0.0, sizeof(float)*512*nffts);
      for (i=0; i<nffts; i++) {           // for each half symbol
         // i = index of FFT
         // prepare a FFT input (of size 512 in-phase, quadrature samples)
         for(j=0; j<512; j++ ) {
            // FFT inputs are stored sequentially, in blocks 512 in-phase,
            // quadrature samples (represented as complex number).
            // Each FFT is over two symbols (512 samples). From FFT-to-FFT,
            // there is a half symbol step in the signal (128 samples).
            // The samples are filtered with a sine filter.
            k=i*128+j;
            fftin[j][0]=idat[k] * w[j];
            fftin[j][1]=qdat[k] * w[j];
         }
         // Compute one-dimensional DFT using PLAN3
         // results are stored in-order in array fftout,
         // with the zero-frequency (DC) component in fftout[0].
         // The data is an array of type fftw_complex, which is by default a
         // double[2] composed of the real (fftout[i][0]) and imaginary
         // (fftout[i][1]) parts of a complex number.
         // Use the standard “in-order” output ordering—the k-th output
         // corresponds to the frequency k/n (or k/T, where T is your total
         // sampling period). For those who like to think in terms of positive
         // and negative frequencies, this means that the positive frequencies
         // are stored in the first half of the output and the negative
         // frequencies are stored in backwards order in the second half of the
         // output. (The frequency -k/n is the same as the frequency (n-k)/n.)
         fftwf_execute(PLAN3);
         // for each of the 512 frequency bins
         for (j=0; j<512; j++ ) {
            k=j+256; // start with negative frequencies
            if( k>511 )
               k=k-512;
            // Calculate the power (V^2) at each frequency
            // Results are stored in matrix "ps" (each row is a FFT, each
            // column is a frequency).
            ps[j][i]=fftout[k][0]*fftout[k][0]+fftout[k][1]*fftout[k][1];
         }
      }
      // Compute the power for each frequency.
      // Result is stored in array "psavg", one cell for each frequency
      memset(psavg,0.0, sizeof(float)*512);
      for (i=0; i<nffts; i++) {
         for (j=0; j<512; j++) {
            psavg[j]=psavg[j]+ps[j][i];
         }
      }
      // Smooth the average power ith 7-point window and limit
      // spectrum to +/-150 Hz.
      // 411/2 * df = 150 Hz
      int window[7]={1,1,1,1,1,1,1};
      float smspec[411];
      for (i=0; i<411; i++) {           // for each bin
         smspec[i]=0.0;
         for(j=-3; j<=3; j++) {
            // start from FFT bin 256-205=51.
            // 51 bins/512 bins * 187.Hz = 18.7 Hz
            k=256-205+i+j;
            smspec[i]=smspec[i]+window[j+3]*psavg[k];
         }
      }
      // Sort spectrum values, then pick off noise level as a percentile
      float tmpsort[411];


      for (j=0; j<411; j++) {
         tmpsort[j]=smspec[j];
      }
      qsort(tmpsort, 411, sizeof(float), floatcomp);

      // Noise level of spectrum is estimated as 123/411= 30'th percentile
      float noise_level = tmpsort[122];

      /* Renormalize spectrum so that (large) peaks represent an estimate of snr.
       * We know from experience that threshold snr is near -7dB in wspr bandwidth,
       * corresponding to -7-26.3=-33.3dB in 2500 Hz bandwidth.
       * The corresponding threshold is -42.3 dB in 2500 Hz bandwidth for WSPR-15. */

      float min_snr, snr_scaling_factor;
      min_snr = pow(10.0,-7.0/10.0);           //this is min snr in wspr bw
      if( wspr_type == 2 ) {
         snr_scaling_factor=26.3;
      } else {
         snr_scaling_factor=35.3;
      }
      for (j=0; j<411; j++) {
         smspec[j]=smspec[j]/noise_level - 1.0;
         if( smspec[j] < min_snr) smspec[j]=0.1*min_snr;
         continue;
      }

      // Find all local maxima in smoothed spectrum.
      for (i=0; i<200; i++) {
         freq0[i]=0.0;
         snr0[i]=0.0;
         drift0[i]=0.0;
         shift0[i]=0;
         sync0[i]=0.0;
      }

      int npk=0;
      unsigned char candidate;
      // THE IF PART NOT USED
      if( more_candidates ) {
         for(j=0; j<411; j=j+2) {
            candidate = (smspec[j]>min_snr) && (npk<200);
            if ( candidate ) {
               freq0[npk]=(j-205)*df;
               snr0[npk]=10*log10(smspec[j])-snr_scaling_factor;
               npk++;
            }
         }
      } else {
         for(j=1; j<410; j++) {
            candidate = (smspec[j]>smspec[j-1]) &&
                        (smspec[j]>smspec[j+1]) &&
                        (npk<200);
            if ( candidate ) {
               freq0[npk]=(j-205)*df; // start from -150 Hz (205*0.73)
               snr0[npk]=10*log10(smspec[j])-snr_scaling_factor;
               npk++;
            }
         }
      }


      // fmin, fmax AND following "for" NOT USED? Because dialfreq_error=0
      // Compute corrected fmin, fmax, accounting for dial frequency error
      fmin += dialfreq_error;           // dialfreq_error is in units of Hz
      fmax += dialfreq_error;

      // Don't waste time on signals outside of the range [fmin,fmax].
      i=0;
      for( j=0; j<npk; j++) {
         if( freq0[j] >= fmin && freq0[j] <= fmax ) {
            freq0[i]=freq0[j];
            snr0[i]=snr0[j];
            i++;
         }
      }
      npk=i;

      // bubble sort on snr, bringing freq along for the ride
      int pass;
      float tmp;
      for (pass = 1; pass <= npk - 1; pass++) {
         for (k = 0; k < npk - pass; k++) {
            if (snr0[k] < snr0[k+1]) {
               tmp = snr0[k];
               snr0[k] = snr0[k+1];
               snr0[k+1] = tmp;
               tmp = freq0[k];
               freq0[k] = freq0[k+1];
               freq0[k+1] = tmp;
            }
         }
      }

      t0=clock();

      /* Make coarse estimates of shift (DT), freq, and drift

       * Look for time offsets up to +/- 8 symbols (about +/- 5.4 s) relative
         to nominal start time, which is 2 seconds into the file

       * Calculates shift relative to the beginning of the file

       * Negative shifts mean that signal started before start of file

       * The program prints DT = shift-2 s

       * Shifts that cause sync vector to fall off of either end of the data
         vector are accommodated by "partial decoding", such that missing
         symbols produce a soft-decision symbol value of 128

       * The frequency drift model is linear, deviation of +/- drift/2 over the
         span of 162 symbols, with deviation equal to 0 at the center of the
         signal vector.
       */
      int idrift,ifr,if0,ifd,k0;
      int kindex;
      float smax,ss,pow,p0,p1,p2,p3;
      for(j=0; j<npk; j++) {                            //For each candidate...
         smax=-1e30;
         if0=freq0[j]/df+256;  // map frequency back to coefficient index
         for (ifr=if0-2; ifr<=if0+2; ifr++) { //Freq search
            // search in the interval [-10 22] half-symbols,
            // i.e., [-5 11] symbols
            // for( k0=-10; k0<22; k0++) { //Time search
            for( k0=0; k0<24; k0++) { //Time search
               for (idrift=-maxdrift; idrift<=maxdrift; idrift++) { //Drift search
                  ss=0.0;
                  pow=0.0;
                  for (k=0; k<162; k++) { //Sum over symbols
                     ifd=ifr+((float)k-81.0)/81.0*( (float)idrift )/(2.0*df);
                     kindex=k0+2*k;

                     if( kindex < nffts ) {
                        p0=ps[ifd-4][kindex];
                        p1=ps[ifd-1][kindex];
                        p2=ps[ifd+1][kindex];
                        p3=ps[ifd+4][kindex];

                        p0=sqrt(p0); // data=0 synchro=0
                        p1=sqrt(p1); // data=0 synchro=1
                        p2=sqrt(p2); // data=1 synchro=0
                        p3=sqrt(p3); // data=1 synchro=1

                        ss=ss+(2*pr3[k]-1)*((p1+p3)-(p0+p2));
                        pow=pow+p0+p1+p2+p3;
                     }
                  }
                  sync1=ss/pow;
                  if( sync1 > smax ) { //Save coarse parameters
                     smax=sync1;
                     shift0[j]=128*(k0+1);
                     drift0[j]=idrift;
                     freq0[j]=(ifr-256)*df;
                     sync0[j]=sync1;
                  }
               }
            }
         }
      }
      tcandidates += (float)(clock()-t0)/CLOCKS_PER_SEC;
      /*
         Refine the estimates of freq, hift using sync as a metric.
         Sync is calculated such that it is a float taking values in the range
         [0.0,1.0].

         Function sync_and_demodulate has three modes of operation
         mode is the last argument:

         0 = no frequency or drift search. find best time lag.
         1 = no time lag or drift search. find best frequency.
         2 = no frequency or time lag search. Calculate soft-decision
         symbols using passed frequency and shift.

         NB: best possibility for OpenMP may be here: several worker threads
         could each work on one candidate at a time.
       */
      for (j=0; j<npk; j++) {
         memset(symbols,0,sizeof(char)*nbits*2);
         memset(callsign,0,sizeof(char)*13);
         memset(call_loc_pow,0,sizeof(char)*23);

         f1=freq0[j];
         drift1=drift0[j];
         shift1=shift0[j];
         sync1=sync0[j];


         // coarse-grid lag and freq search, then if sync>minsync1 continue
         fstep=0.0; ifmin=0; ifmax=0;
         lagmin=shift1-128;
         lagmax=shift1+128;
         lagstep=64;
         t0 = clock();
         sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                             lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);
         tsync0 += (float)(clock()-t0)/CLOCKS_PER_SEC;

         fstep=0.25; ifmin=-2; ifmax=2;
         t0 = clock();
         sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                             lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);

         // refine drift estimate
         fstep=0.0; ifmin=0; ifmax=0;
         float driftp,driftm,syncp,syncm;
         driftp=drift1+0.5;
         sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                             lagmin, lagmax, lagstep, &driftp, symfac, &syncp, 1);

         driftm=drift1-0.5;
         sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                             lagmin, lagmax, lagstep, &driftm, symfac, &syncm, 1);

         if(syncp>sync1) {
            drift1=driftp;
            sync1=syncp;
         } else if (syncm>sync1) {
            drift1=driftm;
            sync1=syncm;
         }

         tsync1 += (float)(clock()-t0)/CLOCKS_PER_SEC;

         // fine-grid lag and freq search
         if( sync1 > minsync1 ) {

            lagmin=shift1-32; lagmax=shift1+32; lagstep=16;
            t0 = clock();
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);
            tsync0 += (float)(clock()-t0)/CLOCKS_PER_SEC;

            // fine search over frequency
            fstep=0.05; ifmin=-2; ifmax=2;
            t0 = clock();
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);
            tsync1 += (float)(clock()-t0)/CLOCKS_PER_SEC;

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

            // Use mode 2 to get soft-decision symbols
            t0 = clock();
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep,
                                &jiggered_shift, lagmin, lagmax, lagstep, &drift1, symfac,
                                &sync1, 2);
            tsync2 += (float)(clock()-t0)/CLOCKS_PER_SEC;

            sq=0.0;
            for(i=0; i<162; i++) {
               y=(float)symbols[i] - 128.0;
               sq += y*y;
            }
            rms=sqrt(sq/162.0);

            if((sync1 > minsync2) && (rms > minrms)) {
               deinterleave(symbols);
               t0 = clock();

               if ( stackdecoder ) {
                  not_decoded = jelinek(&metric, &cycles, decdata, symbols, nbits,
                                        stacksize, stack, mettab,maxcycles);
               } else {
                  not_decoded = fano(&metric,&cycles,&maxnp,decdata,symbols,nbits,
                                     mettab,delta,maxcycles);
               }

               tfano += (float)(clock()-t0)/CLOCKS_PER_SEC;

            }
            idt++;
            if( quickmode ) break;
         }

         // signal is worth a try and FEC (Jelinek or Fano) successful decoding
         if( worth_a_try && !not_decoded ) {
            ndecodes_pass++;
            // copy decoded data into message, modulo 128
            for(i=0; i<11; i++) {
               if( decdata[i]>127 ) {
                  message[i]=decdata[i]-256;
               } else {
                  message[i]=decdata[i];
               }
            }

            // Unpack the decoded message, update the hashtable, apply
            // sanity checks on grid and power, and return
            // call_loc_pow string and also callsign (for de-duping).
            //printf("Unpacking!\n");
            //noprint=unpk_(message,hashtab,call_loc_pow,callsign);

            // subtract even on last pass
            //if( subtraction && (ipass < npasses ) && !noprint ) {
            //   if( get_wspr_channel_symbols(call_loc_pow, hashtab, channel_symbols) ) {
            //      printf("Subtracting!\n");
            //      subtract_signal2(idat, qdat, npoints, f1, shift1, drift1, channel_symbols);
            //   } else {
            //      break;
            //   }
            //}

            // Remove dupes (same callsign and freq within 3 Hz)
            int dupe=0;
            //for (i=0; i<uniques; i++) {
            //   if(!strcmp(callsign,allcalls[i]) &&
            //      (fabs(f1-allfreqs[i]) <3.0)) dupe=1;
            //}
            // if not a dupe within 3 Hz,
            // then increment "uniques" and copy info into"decodes"
            if( (verbose || !dupe) && !noprint) {
               strcpy(allcalls[uniques],callsign);
               allfreqs[uniques]=f1;
               uniques++;

               // Add an extra space at the end of each line so that wspr-x doesn't
               // truncate the power (TNX to DL8FCL!)

               if( wspr_type == 15 ) {
                  freq_print=dialfreq+(1500+112.5+f1/8.0)/1e6;
                  dt_print=shift1*8*dt-1.0;
               } else {
                  freq_print=dialfreq+(1500+f1)/1e6;
                  dt_print=shift1*dt-1.0;
               }

               strcpy(decodes[uniques-1].date,date);
               strcpy(decodes[uniques-1].time,uttime);
               decodes[uniques-1].sync=sync1;
               decodes[uniques-1].snr=snr0[j];
               decodes[uniques-1].dt=dt_print;
               decodes[uniques-1].freq=freq_print;
               //strcpy(decodes[uniques-1].message,call_loc_pow);
               for(i=0; i<7; i++) {
                     decodes[uniques-1].message[i]=message[i];
               }
               decodes[uniques-1].message[7]='\n';
               decodes[uniques-1].drift=drift1;
               decodes[uniques-1].cycles=cycles;
               decodes[uniques-1].jitter=ii;
            }
         }
      }

      // if( ipass == 0 && writec2 ) {
      //   char c2filename[15];
      //   double carrierfreq=dialfreq;
      //   int wsprtype=2;
      //   strcpy(c2filename,"000000_0001.c2");
      //   printf("Writing %s\n",c2filename);
      //   writec2file(c2filename, wsprtype, carrierfreq, idat, qdat);
      //}
   }  // end of main loop

   // decoded messages are "decodes[i].message",
   // i in 0,...,uniques
   // sort the result in order of increasing frequency
   //struct result temp;
   //for (j = 1; j <= uniques - 1; j++) {
   //  for (k = 0; k < uniques - j; k++) {
   //     if (decodes[k].freq > decodes[k+1].freq) {
   //          temp = decodes[k];
   //        decodes[k]=decodes[k+1];;
   //        decodes[k+1] = temp;
   //     }
   //  }
   // }

   // print running time
   gettimeofday(&end, NULL);
   float deltatime = ((end.tv_sec  - start.tv_sec) * 1000000u +
                      end.tv_usec - start.tv_usec) / 1.e6;
   printf("Decoded in : %.2f seconds; uniques: %d\n",deltatime, uniques);

   if ((fp_fftwf_wisdom_file = fopen(wisdom_fname, "w"))) {
      fftwf_export_wisdom_to_file(fp_fftwf_wisdom_file);
      fclose(fp_fftwf_wisdom_file);
   }

   if(writenoise == 999) return -1;       //Silence compiler warning

   // number of decoded messages
   *uniques_ptr = uniques;
   return 0;
}

// End of WSPR Code

// method for decoding thread
void decoder_impl::decode()
{
   float * idat; float * qdat;      // in-phase and quadrature signals
   int i; int j;
   // initialization
   idat=(float *)malloc(sizeof(float)*fl);
   qdat=(float *)malloc(sizeof(float)*fl);
   // size of output buffer (packed, each byte contains eight data bits)
   int nbytes_out = 7;  // (50 bits) last byte uses only 2 most signif. bits
   // decoding results
   struct result decodes[50];
   int uniques; // number of decoded messages

   // main loop
   while (true) {
      j = 0;
      uniques=0;
      // block until work is available
      sem_wait(&semaphore);
      //fprintf(stderr, "complex is: %.9f.%.9f\n",
      //    buffer[0].real(),buffer[0].imag());
      //fprintf(stderr, "<buffer size is: %d\n",(int)buffer.size());
      mutex.lock(); // get exclusive access to buffer
      assert(buffer.size() >= fl);
      // pop nine seconds of data
      for (i = 0; i < 9*fs; i++) {
         // peek element
         idat[j] = buffer.front().real();
         qdat[j] = buffer.front().imag();
         j++;
         // pop element
         buffer.pop_front();
      }
      //fprintf(stderr, "<buffer size is: %d\n",(int)buffer.size());
      // peek 111 seconds of data
      for (i = 0; i < 111*fs; i++) {
         // peek element
         idat[j] = buffer[i].real();
         qdat[j] = buffer[i].imag();
         j++;
      }
      //fprintf(stderr, "*buffer size is: %d\n",(int)buffer.size());
      i = buffer.size();
      mutex.unlock(); // release exclusive access to buffer
      assert(j==fl);
      // WSPR decoding
      if (int e=wsprd(idat, qdat, decodes, &uniques)) {
         // failed to decode
         std::cerr << "WSPR decoding error:" << e <<  "\n";
      } else {
         fprintf(stderr, "Decoded %d messages\n", uniques);
         // push messages on async out port
         out_mutex.lock();
         for (i=0; i<uniques; i++) {
            printf("--- decoder_imp is pushing unpacked message\n");
            // printf("decoder_imp: %-s \n", decodes[i].message);
            // copy a decoded messages into PDU
            std::vector<unsigned char> vec(7);
            for (int j=0; j<7; j++ ) {
              vec[j] = decodes[i].message[j];
            }
            // send the vector
            pmt::pmt_t vecpmt(pmt::make_blob(&vec[0], 7));
            pmt::pmt_t pdu(pmt::cons(pmt::PMT_NIL, vecpmt));
            message_port_pub(out_port, pdu);
         }
         out_mutex.unlock();
      }
      time_t rawtime;
      struct tm *info;
      time( &rawtime );
      info = localtime( &rawtime );
      fprintf(stderr, "Thread completion time: %s\n", asctime(info));
   }
}

int
decoder_impl::work(int noutput_items,
                   gr_vector_const_void_star &input_items,
                   gr_vector_void_star &output_items)
{
   const gr_complex *in = (const gr_complex *) input_items[0];
   //float *out = (float *) output_items[0];
   int npushed;      // number of samples pushed in ring buffer
   time_t rawtime;
   struct tm *info;

   mutex.lock();
   for (npushed = 0; npushed < noutput_items; npushed++) {
      buffer.push_back(in[npushed]);
      //fprintf(stderr, "in complex is: %.9f.%.9f\n",
      //   in[npushed].real(),in[npushed].imag());
      //fprintf(stderr, "buffer complex is: %.9f.%.9f\n",
      //    buffer[npushed].real(),buffer[npushed].imag());
   }
   mutex.unlock();
   count = count + npushed;

   //fprintf(stderr, "count: %d\n",count);
   if (count >= 45000) {      // 120 of samples ready to be process
      // handoff work
      time( &rawtime );
      info = localtime( &rawtime );
      //fprintf(stderr, "complex is: %.9f.%.9f\n",
      //    buffer[0].real(),buffer[0].imag());
      fprintf(stderr, "Handoff time: %s\n", asctime(info));
      // unlock one thread
      sem_post(&semaphore);
      count = count - (9*fs);           // subtract 9 seconds of data
   }

   // Tell runtime system how many output items we produced.
   return npushed;
}

}   /* namespace uwspr */
} /* namespace gr */
