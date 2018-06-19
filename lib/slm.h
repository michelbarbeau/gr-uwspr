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

// Straight Line Model
#ifndef INCLUDED_SLM_H
#define INCLUDED_SLM_H

#include "candidate_t.h"

namespace gr {
 namespace uwspr {
   class SLM
   {
      public:
       // returns frequency drift according to the straight line model
       int slmFrequencyDrift(mode_nonlinear m_nl, int cf, float t);
       // generator of trajectory parameters for the straight line model
       bool slmGenerator(mode_nonlinear *m_nl);
       // initialize the generator
       void slmGeneratorInit();
       // control variable of generator
       int step;
       SLM() {};
       ~SLM() {};
   };
 } // namespace uwspr
} // namespace gr
#endif /* INCLUDED_SLM_H */
