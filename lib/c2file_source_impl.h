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

#ifndef INCLUDED_UWSPR_C2FILE_SOURCE_IMPL_H
#define INCLUDED_UWSPR_C2FILE_SOURCE_IMPL_H

#include <uwspr/c2file_source.h>

namespace gr {
namespace uwspr {

class c2file_source_impl : public c2file_source
{
private:
   bool d_repeat;
   int maxpts; // maximum number of points
   float *idat; float *qdat; // in-phase and quadrature signals
   unsigned sample_idx; // sample index

public:
   c2file_source_impl(const char *filename, bool repeat);
   ~c2file_source_impl();
   // Where all the action really happens
   int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
};

}   // namespace uwspr
} // namespace gr

#endif /* INCLUDED_UWSPR_C2FILE_SOURCE_IMPL_H */
