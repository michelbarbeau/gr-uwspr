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
#include "FDR_impl.h"
#include <stdlib.h>
#include <math.h>

#ifndef DEBUG
#define DEBUG 1 // set debug mode
#endif
#include "debugmacro.h"

namespace gr {
  namespace uwspr {

    FDR::sptr
    FDR::make(int fs, int fl, int spb, int maxdrift, int maxfreqs,
      int halfbandwidth)
    {
      return gnuradio::get_initial_sptr
        (new FDR_impl(fs, fl, spb, maxdrift, maxfreqs, halfbandwidth));
    }

    /* The private constructor
     */
    FDR_impl::FDR_impl(int fs, int fl, int spb, int maxdrift, int maxfreqs,
      int halfbandwidth)
      : gr::block("FDR",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0))
    {
      // create and register the asynchronous input port
      in_port = pmt::mp("in");
      message_port_register_in(in_port);
      // set up a message handler
      set_msg_handler(in_port, boost::bind(&FDR_impl::transform, this, _1));
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
      // init array of candidate frequencies
      candidates = new candidate_t[maxfreqs];
      // half pass bandwidth
      this->halfbandwidth = halfbandwidth;
      // init DFT size (over two symbols)
      size = 2*spb;
      log("DFT size is %d bins", size);
      int maxfreq = (int)((float)fs/2.0);
      log("Max frequency range is plus, minus %d Hertz", maxfreq);
      // validate half band pass bandwidth
      if (halfbandwidth>maxfreq) {
        fprintf(stderr,
          "Half pass bandwidth (%d) must be lower than max freq range (%d)\n",
          halfbandwidth, maxfreq);
          exit(-1);
      }
      log("DFT size (nu) is %d bins", size);
      // init delta frequency
      df = (float)fs/(float)size;
      log("Delta frequency (df) is %9.9f Hz", df);
      tone_1Hz_idx = (int)round(1.0/df);
      tone_3Hz_idx = (int)round(3.0/df);
      log("Tone indices are %d (1 Hz) and %d (3 Hz)",tone_1Hz_idx,tone_3Hz_idx);
      m = size/2;
      log("Frequency index range is -%d,...,-1,0,1,...,%d", m, m-1);
      hpbm = ceil( (int) (float)halfbandwidth/df);
      log("Frequency index range in pass band is -%d,...,-1,0,1,...,%d",
         hpbm, hpbm-1);
      // windowing function (half-sine wave cycle [0 pi]), in "size" steps
      w = new float [size];
      // generate a windowing function (half-sine wave cycle [0 pi]), in "size" steps
      for(int i=0; i<size; i++) {
         w[i]=sin((M_PI/(size-1))*i);
      }
      // Total number of samples ("fl") divided by "size" is the number of two
      // symbol chunks. Four times the number of symbols is the number of half
      // symbol steps (i.e., n).
      n=floor(((float)fl/(float)spb)*2.0)-3;
      log("Number of DFTs (n) is %d", n);
      // power spectrum representation (float ps[n][size])
      ps = new float*[n]; // matrix of "size" rows
      for(int i = 0; i < n; ++i) {
         ps[i] = new float[size]; // row of "ns" elements
      }
      // power spectrum average
      psavg = new float [size];
      // --- preliminaries for one-dimensional DFT of size "size"
      // Do windowed DFTs over two symbols, stepped by half symbols
      // Note: One symbol (baud) is "size/2" samples. Two symbols times
      // is "size" samples.
      // allocate memory for FFT input array
      fftin=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*size);
      // allocate memory for FFT output array
      fftout=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*size);
      // plan allocation (create object containing data needed by FFTW)
      // size = size of transform
      // fftin = pointer to input array
      // fftout = pointer to output array
      // FFTW_FORWARD (-1) = direction of transform (sign of exponent in
      // transform) PATIENCE = FFTW_ESTIMATE
      PLAN3 = fftwf_plan_dft_1d(size,fftin,fftout,FFTW_FORWARD,FFTW_ESTIMATE);
      // Values used to renormalize spectrum so that (large) peaks represent an
      // estimate of SNR. We know from experience that threshold SNR is near
      // -7 dB in WSPR bandwidth (i.e., 6 Hz), corresponding to -7-26.3=-33.3 dB
      // in 2500 Hz bandwidth.
      min_snr = pow(10.0,-7.0/10.0); // min SNR in WSPR bandwidth (6 Hz)
      log("Minimum (6 Hertz) SNR is %9.2f dB", 10*log10(min_snr));
      snr_scaling_factor = 10*log10((float)6/(float)(2*halfbandwidth));
      log("Pass band (%d Hertz) SNR scaling factor is %9.2f dB",
         halfbandwidth*2, snr_scaling_factor);
      // DFT wisdom file
      char wisdom_fname[200];
      strcpy(wisdom_fname,".");
      strncat(wisdom_fname,"/wspr_wisdom.dat",20);
      FILE *fp_fftwf_wisdom_file;
      if ((fp_fftwf_wisdom_file = fopen(wisdom_fname, "r"))) {
         fftwf_import_wisdom_from_file(fp_fftwf_wisdom_file);
         fclose(fp_fftwf_wisdom_file);
      }
    }

    /* virtual destructor.
     */
    FDR_impl::~FDR_impl()
    {
      delete [] w;
      for(int i = 0; i < n; ++i) {
        delete [] ps[i];
      }
      delete [] ps;
      delete [] psavg;
      fftwf_free(fftin);
      fftwf_free(fftout);
      fftwf_destroy_plan(PLAN3);
    }

    int FDR_impl::floatcomp(const void* elem1, const void* elem2)
    {
       if(*(const float*)elem1 < *(const float*)elem2)
          return -1;
       return *(const float*)elem1 > *(const float*)elem2;
    }

    // print a list of frequencies
    int FDR_impl::printfrequencies(int npk, float freq[])
    {
      fprintf(stdout, "Candidate frequencies (%d) are: ", npk);
      for (int i=0; i<npk; i++) {
        fprintf(stdout, " %9.9f", freq[i]);
      }
      fprintf(stdout, "\n");
    }

#include "pr3.h"

    /* Handler of input message
     */
    void FDR_impl::transform(pmt::pmt_t msg)
    {
      int i, j, k;
      // extract input from PDU
      pmt::pmt_t meta(pmt::car(msg)); // CAR is NULL (dictionary)
      // in-phase and quadrature
      pmt::pmt_t vector(pmt::cdr(msg)); // CDR is vector of "fl" pmt:complex
      std::complex<double> c = pmt::to_complex(pmt::vector_ref(vector, 0));
      for (i=0; i<n; i++) { // for each half symbol
         // prepare a DFT input (of "size" 512" in-phase, quadrature samples)
         for(j=0; j<size; j++ ) {
            // from DFT-to-DFT, there is a half symbol step in the signal
            // ("spb"" samples)
            k=i*(spb/2)+j;
            c = pmt::to_complex(pmt::vector_ref(vector, k));
            // samples are filtered with a sine filter.
            fftin[j][0]=c.real() * w[j];
            fftin[j][1]=c.imag() * w[j];
         }
         // Compute one-dimensional DFT using PLAN3
         // results are stored in-order in array fftout,
         // with the zero-frequency (DC) component in fftout[0].
         // The data is an array of type fftw_complex, which is a
         // double[2] composed of the real (fftout[i][0]) and imaginary
         // (fftout[i][1]) parts of a complex number.
         // The output is in the standard “in-order” output ordering. The k-th
         // corresponds to the frequency k/n. This means that the positive
         // frequencies are stored in the first half of the output and negative
         // frequencies are stored in backwards order in the second half of the
         // output. (The frequency -k/n is the same as the frequency (n-k)/n.)
         fftwf_execute(PLAN3);
         // for each of the "size" frequency bins
         for (j=0; j<size; j++) {
            k=j+spb; // start with negative frequencies
            if(k>size-1) k=k-size;
            // Calculate the power (V^2) at each frequency
            // Results are stored in matrix "ps" (each row is a half symbol,
            // each column is  a half symbol)
            ps[i][j]=fftout[k][0]*fftout[k][0]+fftout[k][1]*fftout[k][1];
         }
      }
      // Compute the power for each frequency.
      // Result is stored in array "psavg", one cell for each frequency
      for (j=0; j<size; j++) {
         psavg[j] = 0;
         // sum of all elements on column "j"
         for (i=0; i<n; i++) {
            psavg[j] = psavg[j]+ps[i][j];
         }
      }
      // Smooth the average power, solely within the passs band
      int finpb = 2*hpbm; // number of frequencies in pass band
      log("Number of frequencies in pass band is %d", finpb);
      float smspec[finpb];
      for (i=0; i<finpb; i++) { // for each bin
         smspec[i]=0.0;
         for(j=-3; j<=3; j++) {
            // start from DFT bin "m-hpm"
            k=m-hpbm+i+j;
            smspec[i]=smspec[i]+psavg[k];
         }
      }
      // Sort spectrum values, then pick off noise level as a percentile
      float tmpsort[finpb];
      for (j=0; j<finpb; j++) {
         tmpsort[j]=smspec[j];
      }
      qsort(tmpsort, finpb, sizeof(float), floatcomp);
      // Noise level of spectrum is estimated as 30'th percentile
      int noiseidx = (int)floor(0.3*(float)finpb);
      log("Noise index is %d", noiseidx);
      float noise_level = tmpsort[ noiseidx ];
      log("Noise level is %9.9f", noise_level);
      for (j=0; j<finpb; j++) {
         // SNR calculation in linear form
         smspec[j]=smspec[j]/noise_level - 1.0;
         if(smspec[j] < min_snr) smspec[j]=0.1*min_snr;
      }
      // Find all local maxima in smoothed spectrum.
      int npk=0; // count number of candidate frequencies
      for(j=1; j<(finpb-1); j++) {
            // energy peak?
            if ( (smspec[j]>smspec[j-1]) &&
                 (smspec[j]>smspec[j+1]) &&
                 (npk<maxfreqs ) ) {
               // start from -"half pass bandwidth" Hz
               candidates[npk].freq=(j-hpbm)*df;
               //candidates[npk].snr=10*log10(smspec[j])+snr_scaling_factor;
               // return 6 Hz SNR in dB form
               candidates[npk].snr=10*log10(smspec[j]);
               npk++;
            }
      }
      // bubble sort on SNR
      int pass;
      candidate_t tmp;
      for (pass = 1; pass <= npk - 1; pass++) {
         for (k = 0; k < npk - pass; k++) {
            if (candidates[k].snr < candidates[k+1].snr) {
               tmp = candidates[k];
               candidates[k] = candidates[k+1];
               candidates[k+1] = tmp;
            }
         }
      }
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
      float sync1;
      float smax,ss,pow,p[4];
      for(j=0; j<npk; j++) { // for each candidate...
         smax=-1e30;
         if0=candidates[j].freq/df+m; // frequency to index
         for (ifr=if0-2; ifr<=if0+2; ifr++) { //Freq search
            // search in the interval [1 26] half-symbols,
            for( k0=0; k0<26; k0++) { //Time search
               for (idrift=-maxdrift; idrift<=maxdrift; idrift++) {                          //Drift search
                  ss=0.0;
                  pow=0.0;
                  for (k=0; k<162; k++) { // sum of power over symbols
                     ifd=ifr+((float)k-81.0)/81.0*( (float)idrift )/(2.0*df);
                     kindex=k0+2*k;
                     p[0]=sqrt(ps[kindex][ifd-tone_3Hz_idx]);
                     p[1]=sqrt(ps[kindex][ifd-tone_1Hz_idx]);
                     p[2]=sqrt(ps[kindex][ifd+tone_1Hz_idx]);
                     p[3]=sqrt(ps[kindex][ifd+tone_3Hz_idx]);
                     ss=ss+(2*pr3[k]-1)*((p[1]+p[3])-(p[0]+p[2]));
                     pow=pow+p[0]+p[1]+p[2]+p[3];
                  }
                  sync1=ss/pow;
                  if(sync1 > smax) {                               //Save coarse parameters
                     smax=sync1;
                     candidates[j].shift=128*(k0+1);
                     candidates[j].drift=idrift;
                     candidates[j].freq=(ifr-m)*df; // index to frequency
                     candidates[j].sync=sync1;
                     // a linear case
                     candidates[j].m_type = linear;
                     candidates[j].m_linear.drift = idrift;
                  }
               }
            }
         }
      }
#if DEBUG
      // printfrequencies(npk, freq0);
#endif
      // 2nd tuple component is number of candidate frequencies
      pmt::pmt_t n = pmt::from_long(npk);
      // 3rd tuple component is frequencies
      pmt::pmt_t frequencies = pmt::make_vector(npk, pmt::from_double(0.0));
      // 4th tuple component is SNRs
      pmt::pmt_t SNRs = pmt::make_vector(npk, pmt::from_double(0.0));
      // 5th tuple componet is drifts
      pmt::pmt_t DRIFTs = pmt::make_vector(npk, pmt::from_double(0.0));
      // 6th tuple componet is syncs
      pmt::pmt_t SYNCs = pmt::make_vector(npk, pmt::from_double(0.0));
      // 7th tuple componet is syncs
      pmt::pmt_t SHIFTs = pmt::make_vector(npk, pmt::from_long(0));
      for (i = 0; i < npk; i++) {
          pmt::vector_set(frequencies,i,pmt::from_double(candidates[i].freq));
          pmt::vector_set(SNRs,i,pmt::from_double(candidates[i].snr));
          pmt::vector_set(DRIFTs,i,pmt::from_double(candidates[i].drift));
          pmt::vector_set(SYNCs,i,pmt::from_double(candidates[i].sync));
          pmt::vector_set(SHIFTs,i,pmt::from_long(candidates[i].shift));
      }
      // PDU payload is a tuple: 1st compunent is vector of samples (in-phase & quadrature),
      pmt::pmt_t tuple =
         pmt::make_tuple(vector,n,frequencies,SNRs,DRIFTs,SYNCs,SHIFTs);
      // store into PDU: CAR is NULL dictionary, CDR is a tuple
      pmt::pmt_t pdu(pmt::cons(pmt::PMT_NIL, tuple));
      // post the PDU
      message_port_pub(out_port, pdu);
    }

  } /* namespace uwspr */
} /* namespace gr */
