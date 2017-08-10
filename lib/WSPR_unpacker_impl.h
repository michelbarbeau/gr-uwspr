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

#ifndef INCLUDED_UWSPR_WSPR_UNPACKER_IMPL_H
#define INCLUDED_UWSPR_WSPR_UNPACKER_IMPL_H

#include <uwspr/WSPR_unpacker.h>
#include <stdio.h>
#include <stdlib.h>
#include "WSPR_unpacker_impl.h"
#include "helpers.h"

namespace gr {
  namespace uwspr {

    class WSPR_unpacker_impl : public WSPR_unpacker, public helpers
    {
     private:
      pmt::pmt_t in_port;
      pmt::pmt_t out_port;
      void unpack(pmt::pmt_t msg);
      FILE *fp_fftwf_wisdom_file, *fall_wspr, *fwsprd, *fhash, *ftimer;
      char *callsign, *call_loc_pow;
      char *hashtab;
      int usehashtable;
     public:
      WSPR_unpacker_impl();
      ~WSPR_unpacker_impl();

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace uwspr
} // namespace gr

#endif /* INCLUDED_UWSPR_WSPR_UNPACKER_IMPL_H */
