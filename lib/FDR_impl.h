/* -*- c++ -*- */
/*
 * Copyright 2019 Michel Barbeau, Carleton University.
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
#include "candidate_t.h"
#include "slm.h"

namespace gr {
  namespace uwspr {
    class FDR_impl : public FDR, public SLM
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
       // carrier frequency (used for Doppler shift estimation)
       int carrierfrequency;
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
       // array of candidate frequencies
       candidate_t * candidates;
       // used for spectrum normalization
       float min_snr, snr_scaling_factor;
       // print a list of frequencies
       static int printfrequencies(int npk, float freq[]);
       void powersum(
         int k0, int k, int ifd, float *ss, float *pow);
       // Nonlinear vs linear power ratio required
       float threshold;
     public:
      FDR_impl(int fs, int fl, int spb, int maxdrift, int maxfreqs,
        int halfbandwidth, int cf, int threshold);
      ~FDR_impl();
    };
  } // namespace uwspr
} // namespace gr

#endif /* INCLUDED_UWSPR_FDR_IMPL_H */
