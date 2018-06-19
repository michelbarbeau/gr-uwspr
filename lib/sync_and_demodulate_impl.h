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

#ifndef INCLUDED_UWSPR_SYNC_AND_DEMODULATE_IMPL_H
#define INCLUDED_UWSPR_SYNC_AND_DEMODULATE_IMPL_H

#include <uwspr/sync_and_demodulate.h>
#include <fstream>
#include <iostream>
#include "candidate_t.h"
#include "slm.h"

using namespace std;

extern unsigned char Partab[];

namespace gr {
  namespace uwspr {
    class sync_and_demodulate_impl : public sync_and_demodulate, public SLM
    {
     private:
       pmt::pmt_t in_port;
       pmt::pmt_t out_port;
       int printfrequencies(int npk, float freq[]);
       void demodulate(pmt::pmt_t msg);
       // sampling frequency (samples/second)
       int fs;
       // frame length, in samples
       int fl;
       // number of samples by baud (symbol)
       int spb;
       // maximum frequency drift
       int maxdrift;
       // maximum number of frequencies consider during the search
       int maxfreqs;
       // delta frequency
       float df; // or df=775/256
       int npoints;
       int nffts;
       float tfano,treadwav,tcandidates,tsync0;
       float tsync1,tsync2,ttotal;
       unsigned int nbits; // number of bits in a frame
       unsigned char *symbols, *decdata, *channel_symbols;
       // array of candidate frequencies
       candidate_t * candidates;
       // carrier frequency (used for Doppler shift estimation)
       int carrierfrequency;
       // message log file
       FILE * msglogfile;
       // log a message
       void logmessage(signed char message[]);
       // for time stamps
       time_t start;
       struct tm *info;
       void printtime();
       // frame counter
       unsigned framecount;
       // Fano FEC structure
       struct node {
          unsigned long encstate; /* Encoder state of next node */
          long gamma; /* Cumulative metric to this node */
          int metrics[4]; /* Metrics indexed by all possible tx syms */
          int tm[2]; /* Sorted metrics for current hypotheses */
          int i; /* Current branch being tested */
       };
       // metric table
       float bias; //Fano metric bias (used for both Fano and stack algorithms)
       int mettab[2][256];
       // Convolutional FEC (Fano) methods
       int encode(
          unsigned char *symbols,          // Output buffer, 2*nbytes*8
          unsigned char *data,           // Input buffer, nbytes
          unsigned int nbytes);            // Number of bytes in data
       // FANO FEC decoder
       int fano(unsigned int *metric, unsigned int *cycles, unsigned int *maxnp,
          unsigned char *data,unsigned char *symbols, unsigned int nbits,
            int mettab[2][256],int delta,unsigned int maxcycles);
       // synchronize and demodulate
       void sync_and_demodulate(candidate_t candidate,
         float *id, float *qd, long np,
         unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
         int *shift1, int lagmin, int lagmax, int lagstep,
         float *drift1, int symfac, float *sync, int mode);
       void deinterleave(unsigned char *sym);
       static int floatcomp(const void* elem1, const void* elem2);
     public:
      sync_and_demodulate_impl(int fs, int fl, int spb, int maxdrift,
        int maxfreqs, int cf);
      ~sync_and_demodulate_impl();
    };
  } // namespace uwspr
} // namespace gr
#endif /* INCLUDED_UWSPR_SYNC_AND_DEMODULATE_IMPL_H */
