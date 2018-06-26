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

#include <Fano.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef DEBUG
#define DEBUG 0 // set debug mode
#endif
#include "debugmacro.h"

namespace gr {
 namespace uwspr {
   /*
    * The private constructor
    */
   Fano::Fano()
   {
     // metric table initialization
#include "./metric_tables.c"
     bias=0.45;    //Fano metric bias (used for both Fano and stack algorithms)
     for(int i=0; i<256; i++) {
       mettab[0][i]=round( 10*(metric_tables[2][i]-bias) );
       mettab[1][i]=round( 10*(metric_tables[2][255-i]-bias) );
     }
   }

   Fano::~Fano()
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
    * Storing only one symbol per byte uses more space, but it is faster
    * and easier than trying to pack them more compactly.
    */
   int Fano::encode(
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
   int Fano::fano(
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
#if DEBUG
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
  }   /* namespace uwspr */
} /* namespace gr */
