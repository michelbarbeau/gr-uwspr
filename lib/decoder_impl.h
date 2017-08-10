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

#ifndef INCLUDED_UWSPR_DECODER_IMPL_H
#define INCLUDED_UWSPR_DECODER_IMPL_H

#include <uwspr/decoder.h>
#include <boost/circular_buffer.hpp>
#include <semaphore.h>
#include <fftw3.h>
#include "helpers.h"

extern unsigned char Partab[];

namespace gr {
   namespace uwspr {

typedef struct buffer buffer_t;

class decoder_impl : public decoder, public helpers
{
private:
   pmt::pmt_t out_port; // async output port
   boost::mutex out_mutex; // async output port
   void decode(); // method for decoding threads
   // constants
   int fs; // sampling frequency (samples/second)
   int fl; // data length, in samples
   // variables
   boost::circular_buffer<gr_complex> buffer; // circular buffer for in samples
   int count; // number of samples ready to be processed in circular buffer
   boost::mutex mutex; // for circular buffer access
   sem_t semaphore; // controls the pool of threads

// Convolutional FEC (Fano) methods
int encode(
   unsigned char *symbols,          // Output buffer, 2*nbytes*8
   unsigned char *data,           // Input buffer, nbytes
   unsigned int nbytes);            // Number of bytes in data
int fano(unsigned int *metric, unsigned int *cycles, unsigned int *maxnp,
         unsigned char *data,unsigned char *symbols, unsigned int nbits,
         int mettab[2][256],int delta,unsigned int maxcycles);
// WSPR structures
struct node {
   unsigned long encstate;          /* Encoder state of next node */
   long gamma;            /* Cumulative metric to this node */
   int metrics[4];            /* Metrics indexed by all possible tx syms */
   int tm[2];             /* Sorted metrics for current hypotheses */
   int i;               /* Current branch being tested */
};
struct snode {
   uint64_t encstate;        // Encoder state
   int gamma;                    // Cumulative metric to this node
   unsigned int depth;                      // depth of this node
   unsigned int jpointer;
};
// decoding result
struct result { char date[7]; char time[5]; float sync; float snr;
                float dt; double freq; char message[23]; float drift;
                unsigned int cycles; int jitter; };
// WSPR related methods
void sync_and_demodulate(float *id, float *qd, long np,
                         unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
                         int *shift1, int lagmin, int lagmax, int lagstep,
                         float *drift1, int symfac, float *sync, int mode);
void subtract_signal2(float *id, float *qd, long np,
                      float f0, int shift0, float drift0, unsigned char* channel_symbols);
unsigned long writec2file(char *c2filename, int trmin,
                          double freq, float *idat, float *qdat);
int jelinek(
   unsigned int *metric,        /* Final path metric (returned value) */
   unsigned int *cycles,        /* Cycle count (returned value) */
   unsigned char *data,         /* Decoded output data */
   unsigned char *symbols,        /* Raw deinterleaved input symbols */
   unsigned int nbits,        /* Number of output bits */
   unsigned int stacksize,
   struct snode *stack,
   int mettab[2][256],        /* Metric table, [sent sym][rx symbol] */
   unsigned int maxcycles);       /* Decoding timeout in cycles per bit */
static int floatcomp(const void* elem1, const void* elem2);
void deinterleave(unsigned char *sym);

char get_locator_character_code(char ch);
char get_callsign_character_code(char ch);
long unsigned int pack_grid4_power(char *grid4, int power);
long unsigned int pack_call(char *callsign);
void pack_prefix(char *callsign, int32_t *n, int32_t *m, int32_t *nadd );
void interleave(unsigned char *sym);
int get_wspr_channel_symbols(char* message, char* hashtab,
                             unsigned char* symbols);
int wsprd(float * idat, float * qdat, struct result decodes[], int * uniques);
//--- WSPR related data members
struct snode *stack;    // stack
// data files
char *data_dir=NULL;
unsigned long nr;
float tfano,treadwav,tcandidates,tsync0;
float tsync1,tsync2,ttotal;
FILE *fp_fftwf_wisdom_file;
char wisdom_fname[200];
// stack decoder_impl
int stacksize;
int stackdecoder;
// metric table
float bias;    //Fano metric bias (used for both Fano and stack algorithms)
int mettab[2][256];
// date and time
char date[7], uttime[5];
// FFT
unsigned int npoints;
fftwf_complex *fftin, *fftout;
int nffts;
// Possible PATIENCE options: FFTW_ESTIMATE, FFTW_ESTIMATE_PATIENT,
// FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
#define PATIENCE FFTW_ESTIMATE
fftwf_plan PLAN3;
// windowing function (half-sine wave cycle [0 pi]), in 512 steps
float w[512];
// Maximum (+/-) drift
int maxdrift;
public:
decoder_impl(int maxdrift);
~decoder_impl();
// function implementing the main loop
int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
};
}   // namespace uwspr
} // namespace gr

#endif /* INCLUDED_UWSPR_DECODER_IMPL_H */
