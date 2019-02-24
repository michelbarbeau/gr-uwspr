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

#ifndef INCLUDED_CANDIDATE_T_H
#define INCLUDED_CANDIDATE_T_H

namespace gr {
  namespace uwspr {
    // candidate frequency structures
    typedef struct {
      float drift;
    } mode_linear;
    typedef struct {
      // velocity
      double V1, V2;
      // start coordinates
      int p1, p2;
    } mode_nonlinear;
    enum Modes { linear, nonlinear };
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
  } // namespace uwspr
} // namespace gr

#endif /* INCLUDED_CANDIDATE_T_H */
