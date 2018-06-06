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

#ifndef INCLUDED_UWSPR_FDR_IMPL_H
#define INCLUDED_UWSPR_FDR_IMPL_H

#include <uwspr/FDR.h>
#include <fftw3.h>

namespace gr {
  namespace uwspr {

    // candidate frequency structures
    typedef struct {
      float drift;
    } mode_linear;
    typedef struct {
      double va; // velocity
      double ma; // slope
      int x, y; // start coordinates
    } mode_nonlinear;
    enum Modes { linear, nonlinear};
    typedef struct {
       float freq;
       float snr;
       float drift;
       float sync;
       int shift;
       Modes m_type; // mode type (1: linear, 2: nonlinear)
       union {
         // linear frequency drift model
         mode_linear m_linear;
         // straight-line nonlinear frequency drift model
         mode_nonlinear m_nonlinear;
       };
    }  candidate_t;

    class FDR_impl : public FDR
    {
     private:
       // asynchronous input port
       pmt::pmt_t in_port;
       // asynchronous input port handler
       void transform(pmt::pmt_t msg);
       // asynchronous output port
       pmt::pmt_t out_port;
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
       // number of points (samples) used for the DFT
       int npoints;
       // DFT size
       int size;
       // max (baseband) frequency
       int maxfreq;
       // frequency index upper bound
       int m;
       // frequency index upper bound in band pass
       int hpbm;
       // half pass Bandwidth
       int halfbandwidth;
       // pointers to the DFTs
       fftwf_complex *fftin, *fftout;
       // number of DFTs
       int n;
       // DFT plan
       fftwf_plan PLAN3;
       // windowing function (half-sine wave cycle [0 pi]), in 512 steps
       float *w;
       // power spectrum representation
       float **ps;
       // power spectrum average
       float *psavg;
       // min and max frequency range
       float fmin, fmax;
       // float comparison function
       static int floatcomp(const void* elem1, const void* elem2);
       // delta frequency
       float df;
       // indices of modulation tones
       int tone_1Hz_idx, tone_3Hz_idx;
       // array of candidate printfrequencies
       candidate_t * candidates;
       // used for spectrum normalization
       float min_snr, snr_scaling_factor;
       // print a list of frequencies
       static int printfrequencies(int npk, float freq[]);
     public:
      FDR_impl(int fs, int fl, int spb, int maxdrift, int maxfreqs,
        int halfbandwidth);
      ~FDR_impl();
    };
  } // namespace uwspr
} // namespace gr

#endif /* INCLUDED_UWSPR_FDR_IMPL_H */
