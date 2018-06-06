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
#include <iostream>
#include <assert.h>
#include <stdbool.h>
#include <volk/volk.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>

#ifndef DEBUG
#define DEBUG 1 // set debug mode
#endif
#include "debugmacro.h"

namespace gr {
namespace uwspr {

sync_and_demodulate::sptr
sync_and_demodulate::make(int fs, int fl, int spb, int maxdrift)
{
        return gnuradio::get_initial_sptr
                       (new sync_and_demodulate_impl(fs, fl, spb, maxdrift));
}

/*
 * The private constructor
 */
sync_and_demodulate_impl::sync_and_demodulate_impl(
        int fs, int fl, int spb, int maxdrift)
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
}


/*
 * Virtual destructor.
 */
sync_and_demodulate_impl::~sync_and_demodulate_impl()
{
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
 *
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
 *
 * Decode packet with the Fano algorithm.
 * Return 0 on success, -1 on timeout
 *
 * Decode packet with the Fano algorithm.
 * Return 0 on success, -1 on timeout
 */
int sync_and_demodulate_impl::fano(
        unsigned int  *metric,     // Final path metric (returned value)
        unsigned int  *cycles,     // Cycle count (returned value)
        unsigned int  *maxnp,     // Progress before timeout (returned value)
        unsigned char *data,     // Decoded output data
        unsigned char *symbols,     // Raw deinterleaved input symbols
        unsigned int nbits,     // Number of output bits
        int mettab[2][256],     // Metric table, [sent sym][rx symbol]
        int delta,     // Threshold adjust parameter
        unsigned int maxcycles)     // Decoding timeout in cycles per bit
{
        struct node *nodes;     // First node
        struct node *np;           // Current node
        struct node *lastnode;     // Last node
        struct node *tail;     // First node of tail
        int t;       // Threshold
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
        ENCODE(lsym,np->encstate);     // 0-branch (LSB is 0)
        m0 = np->metrics[lsym];

        /* Now do the 1-branch. To save another ENCODE call here and
         * inside the loop, we assume that both polynomials are odd,
         * providing complementary pairs of branch symbols.

         * This code should be modified if a systematic code were used.
         */

        m1 = np->metrics[3^lsym];
        if(m0 > m1) {
                np->tm[0] = m0;                     // 0-branch has better metric
                np->tm[1] = m1;
        } else {
                np->tm[0] = m1;                     // 1-branch is better
                np->tm[1] = m0;
                np->encstate++;                 // Set low bit
        }
        np->i = 0;                                // Start with best branch
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
                        if(np->gamma < t + delta) {     // Node is acceptable
                                /* First time we've visited this node;
                                 * Tighten threshold.
                                 *
                                 *   t += delta * ((ngamma - t)/delta);
                                 * but the multiply and divide are slower.
                                 */
                                while(ngamma >= t + delta) t += delta;
                        }
                        np[1].gamma = ngamma;       // Move forward
                        np[1].encstate = np->encstate << 1;
                        if( ++np == (lastnode+1) ) {
                                break;        // Done!
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
                                        np->tm[0] = m0;     // 0-branch is better
                                        np->tm[1] = m1;
                                } else {
                                        np->tm[0] = m1;     // 1-branch is better
                                        np->tm[1] = m0;
                                        np->encstate++;     // Set low bit
                                }
                        }
                        np->i = 0;            // Start with best branch
                        continue;
                }
                // Threshold violated, can't go forward
                for(;; ) {                          // Look backward
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
                                np->i++;     // Search next best branch
                                np->encstate ^= 1;
                                break;
                        }                     // else keep looking back
                }
        }
        *metric =  np->gamma;                 // Return the final path metric

        // Copy decoded data to user's buffer
        nbits >>= 3;
        np = &nodes[7];
        while(nbits-- != 0) {
                *data++ = np->encstate;
                np += 8;
        }
        *cycles = i+1;

        free(nodes);
        if(i >= maxcycles) return -1;         // Decoder timed out
        return 0;                   // Successful completion
}

#include "pr3.h"

void sync_and_demodulate_impl::sync_and_demodulate(
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
                                fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0; // drift search
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

// print a list of frequencies
int sync_and_demodulate_impl::printfrequencies(int npk, float freq[])
{
        fprintf(stdout, "Candidate frequencies (%d) are: ", npk);
        for (int i=0; i<npk; i++) {
                fprintf(stdout, " %9.9f", freq[i]);
        }
        fprintf(stdout, "\n");
}

void sync_and_demodulate_impl::demodulate(pmt::pmt_t msg)
{
        /* Local variables
         */
        int i,j,k;
        signed char message[]={-9,13,-35,123,57,-39,64,0,0,0,0};
        int delta,maxpts=65536;
        int shift1, lagmin, lagmax, lagstep, ifmin, ifmax, worth_a_try, not_decoded;
        unsigned int metric, cycles, maxnp;
        float df=375.0/256.0/2; // or df=775/256
        float dt=1.0/375.0, dt_print;
        double dialfreq_cmdline=0.0, dialfreq, freq_print;
        double dialfreq_error=0.0;
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
        log("new PDU");
        pmt::pmt_t meta(pmt::car(msg)); // CAR is NULL (dictionary)
        pmt::pmt_t tuple(pmt::cdr(msg)); // CDR is a tuple
        // 1st tuple component is signal time domain representation
        pmt::pmt_t vector(pmt::tuple_ref(tuple,0));
        std::complex<double> c;
        c = pmt::to_complex(pmt::vector_ref(vector,0));
        float idat[fl]; float qdat[fl];
        for(j=0; j<fl; j++) {
                c = pmt::to_complex(pmt::vector_ref(vector,j));
                idat[j] = c.real();
                qdat[j] = c.imag();
        }
        // 2nd tuple component is number of candidate frequencies
        pmt::pmt_t n(pmt::tuple_ref(tuple,1));
        int npk = pmt::to_long(n);
        log("npk = %d",npk);
        // 3rd tuple component is frequencies
        pmt::pmt_t frequencies(pmt::tuple_ref(tuple,2));
        float freq0[200];
        // 4th tuple component is SNRs
        pmt::pmt_t SNRs(pmt::tuple_ref(tuple,3));
        float snr0[200];
        // 5th tuple componet is drifts
        pmt::pmt_t DRIFTs(pmt::tuple_ref(tuple,4));
        float drift0[200];
        // 6th tuple componet is syncs
        pmt::pmt_t SYNCs(pmt::tuple_ref(tuple,5));
        // 7th tuple componet is syncs
        pmt::pmt_t SHIFTs(pmt::tuple_ref(tuple,6));
        int shift0[200];
        float sync0[200];
        for(j=0; j<npk; j++) {
                freq0[j] = pmt::to_double(pmt::vector_ref(frequencies,j));
                snr0[j] = pmt::to_double(pmt::vector_ref(SNRs,j));
                drift0[j] = pmt::to_double(pmt::vector_ref(DRIFTs,j));
                sync0[j] = pmt::to_double(pmt::vector_ref(SYNCs,j));
                shift0[j] = pmt::to_long(pmt::vector_ref(SHIFTs,j));
        }
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
                f1=freq0[j];
                drift1=drift0[j];
                shift1=shift0[j];
                sync1=sync0[j];
                // coarse-grid lag and freq search, then if sync>minsync1 continue
                fstep=0.0; ifmin=0; ifmax=0;
                lagmin=shift1-128;
                lagmax=shift1+128;
                lagstep=64;
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                    lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);
                fstep=0.25; ifmin=-2; ifmax=2;
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

                //tsync1 += (float)(clock()-t0)/CLOCKS_PER_SEC;

                // fine-grid lag and freq search
                if( sync1 > minsync1 ) {

                        lagmin=shift1-32; lagmax=shift1+32; lagstep=16;
                        //t0 = clock();
                        sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                            lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0);
                        //tsync0 += (float)(clock()-t0)/CLOCKS_PER_SEC;

                        // fine search over frequency
                        fstep=0.05; ifmin=-2; ifmax=2;
                        //t0 = clock();
                        sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                            lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1);
                        //tsync1 += (float)(clock()-t0)/CLOCKS_PER_SEC;

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
                        sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep,
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
                        log("*** pushing packed message ***");
                        log("Baseband freq is %9.2f Hz", freq0[j]);
                        log("(6 Hz) SNR is %9.2f dB", snr0[j]);
                        log("Drift is %9.2f Hz", drift0[j]);
                        // push the packed message
                        pmt::pmt_t payload(pmt::make_blob(&message[0], 7));
                        pmt::pmt_t pdu(pmt::cons(pmt::PMT_NIL, payload));
                        message_port_pub(out_port, pdu);
                }
        }
        return;
}
}   /* namespace uwspr */
} /* namespace gr */
