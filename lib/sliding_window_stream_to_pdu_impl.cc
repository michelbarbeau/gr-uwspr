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
#include "sliding_window_stream_to_pdu_impl.h"
#include <time.h>

#ifndef DEBUG
#define DEBUG 0 // set debug mode
#endif
#include "debugmacro.h"

namespace gr {
  namespace uwspr {

    sliding_window_stream_to_pdu::sptr
    sliding_window_stream_to_pdu::make(int fs, int fl, int shift, int C)
    {
      return gnuradio::get_initial_sptr
        (new sliding_window_stream_to_pdu_impl(fs, fl, shift, C));
    }

    /*
     * The private constructor
     */
    sliding_window_stream_to_pdu_impl::sliding_window_stream_to_pdu_impl(
      int fs, int fl, int shift, int C)
      : gr::sync_block("sliding_window_stream_to_pdu",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(0, 0, 0))
    {
      // create and register the asynchronous output port
      out_port = pmt::mp("out");
      message_port_register_out(out_port);
      // init sampling frequency
      this->fs = fs;
      // init frame length, in samples
      this->fl = fl;
      // buffer window shift from PUD to PDU (seconds)
      this->shift = shift;
      // init number of items in ring buffer
      count = 0;
      // init capacity of ring buffer
      buffer.set_capacity(C*fl);
      time(&start);
      struct tm *info = localtime(&start);
      log("Start time: %s", asctime(info));
    }

    /*
     * Our virtual destructor.
     */
    sliding_window_stream_to_pdu_impl::~sliding_window_stream_to_pdu_impl()
    {
    }

    // print current time and ellapsed timed
    void sliding_window_stream_to_pdu_impl::printtime()
    {
      time_t now;
      struct tm *info;
      // handoff work to next block
      time(&now);
      info = localtime(&now);
      log("Handoff time: %s", asctime(info));
      time_t dtime = difftime(now,start);
      int h = (dtime / 3600) % 24;
      int m = (dtime / 60) % 60;
      int s = dtime % 60;
      log("Ellapsed time: %02d:%02d:%02d", h, m, s);
    }

    /*
     * Main loop
     */
    int
    sliding_window_stream_to_pdu_impl::work(int noutput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
       const gr_complex *in = (const gr_complex *) input_items[0];
       int npushed; // number of samples pushed in ring buffer
       int i; // index on PDU samples
       std::complex<double> c(0,0);
       pmt::pmt_t vector = pmt::make_vector(fl, pmt::from_complex(c));
       // store input samples into ring buffer
       for (npushed = 0; npushed < noutput_items; npushed++) {
          buffer.push_back(in[npushed]);
       }
       count = count + npushed;
       // fprintf(stderr, "count: %d\n",count);
       if (count >= fl) { // 120 of samples ready to be process
          // handoff work to next block
          printtime();
          // pop "shift*fs" samples and store into PDU
          for (i = 0; i < shift*fs; i++) {
             // peek
              pmt::vector_set(vector,i,
                pmt::make_rectangular(
                  buffer.front().real(),buffer.front().imag()));
             // pop
             buffer.pop_front();
          }
          // peek "fl-(shift*fs)" samples and store into PDU
          for (i = 0; i < (fl-(shift*fs)); i++) {
             pmt::vector_set(vector,shift*fs+i,
               pmt::make_rectangular(buffer[i].real(),buffer[i].imag()));
          }
          // store into PDU: CAR is NULL dictionary, CDR vector of samples
          pmt::pmt_t pdu(pmt::cons(pmt::PMT_NIL, vector));
          // post the PDU
          message_port_pub(out_port, pdu);
          count = count - (shift*fs); // subtract 9 seconds of data
       }
       // Tell runtime system how many output items we produced.
       return npushed;
    }

  } /* namespace uwspr */
} /* namespace gr */
