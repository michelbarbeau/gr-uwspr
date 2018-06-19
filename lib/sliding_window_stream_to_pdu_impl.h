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

#ifndef INCLUDED_UWSPR_SLIDING_WINDOW_STREAM_TO_PDU_IMPL_H
#define INCLUDED_UWSPR_SLIDING_WINDOW_STREAM_TO_PDU_IMPL_H

#include <uwspr/sliding_window_stream_to_pdu.h>
#include <boost/circular_buffer.hpp>

namespace gr {
  namespace uwspr {
    class sliding_window_stream_to_pdu_impl :
     public sliding_window_stream_to_pdu
    {
     private:
      // asynchronous output port
      pmt::pmt_t out_port;
      // sampling frequency (samples/second)
      int fs;
      // frame length, in samples
      int fl;
      // buffer window shift from PUD to PDU (seconds)
      int shift;
      // circular buffer for input samples
      boost::circular_buffer< std::complex<double> > buffer;
      // number of samples ready to be processed in circular buffer
      int count;
      // for time stamps
      time_t start;
      struct tm *info;
     public:
      // constructor
      sliding_window_stream_to_pdu_impl(int fs, int fl, int shift, int C);
      // destructor
      ~sliding_window_stream_to_pdu_impl();
      // print current time and ellapsed timed
      void printtime();
      // main loop
      int work(int noutput_items,
       gr_vector_const_void_star &input_items,
       gr_vector_void_star &output_items);
    };
  } // namespace uwspr
} // namespace gr

#endif /* INCLUDED_UWSPR_SLIDING_WINDOW_STREAM_TO_PDU_IMPL_H */
